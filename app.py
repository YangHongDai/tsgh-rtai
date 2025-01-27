import os
import re
import json
import logging
import requests
from datetime import datetime
import time
from flask import Flask, request
from diskcache import Cache
from linebot.v3 import WebhookHandler
from linebot.v3.messaging import (
    Configuration, 
    ApiClient, 
    MessagingApi, 
    ReplyMessageRequest, 
    TextMessage, FlexMessage, 
    MessagingApiBlob, 
    RichMenuSize, 
    RichMenuRequest, 
    RichMenuArea, 
    RichMenuBounds, 
    MessageAction,
    QuickReply,
    QuickReplyItem,
    PostbackAction,
    ReplyMessageRequest,
    URIAction)
from linebot.v3.webhooks import MessageEvent, TextMessageContent, PostbackEvent
from dotenv import load_dotenv


# ------------------------- 初始化配置 -------------------------
load_dotenv()
app = Flask(__name__)
app.config.update(JSON_AS_ASCII=False)

# 配置參數
DEEPSEEK_API_URL = "https://api.deepseek.com/v1/chat/completions"
CACHE_TTL = 3600  # 1小時快取
MAX_RETRIES = 3    # API呼叫重試次數

# 初始化組件
LINE_CHANNEL_TOKEN = 'cd1N2IYrGKTouMLPWRbgmDUl2DjyHEhDucB/9BGXaKUEWHeiSdc+iKY4v6fMUhZm1cV+bSCJm5uy+H2ZvkJwNiOmixiEqyh5DKbUAsGZFr67xn1VwDwiPP0uGt7dUAJiKhmmxdxyWEa+Fc986K2qgQdB04t89/1O/w1cDnyilFU='
cache = Cache("response_cache")
configuration = Configuration(access_token=os.getenv("LINE_CHANNEL_TOKEN"))
handler = WebhookHandler(os.getenv("LINE_CHANNEL_SECRET"))

# 配置日誌
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ------------------------- 醫療安全檢查模組 -------------------------
class MedicalSafety:
    """醫療安全過濾器"""
    
    def __init__(self):
        self.sensitive_keywords = {
            'emergency': ['心臟病發作', '中風', '大出血', '呼吸困難', '意識喪失'],
            'dangerous': ['自殺', '自殘', '謀殺', '下毒'],
            'sensitive': ['性病', 'HIV', '墮胎', '毒品']
        }
        self.logger = logging.getLogger('MedicalSafety')

    def check_input(self, text):
        """輸入內容安全檢查"""
        text_lower = text.lower()
        
        # 緊急情況檢測
        for keyword in self.sensitive_keywords['emergency']:
            if keyword in text_lower:
                self.logger.warning(f"檢測到緊急情況關鍵字: {keyword}")
                return {
                    "safe": False,
                    "message": f"【三軍總醫院衛教機器人阿泰提醒】檢測到緊急醫療情況關鍵字「{keyword}」，請立即撥打119或前往最近急診室！"
                }

        # 危險行為檢測
        for keyword in self.sensitive_keywords['dangerous']:
            if keyword in text_lower:
                self.logger.warning(f"檢測到危險行為關鍵字: {keyword}")
                return {
                    "safe": False,
                    "message": "【三軍總醫院衛教機器人阿泰提醒】檢測到潛在危險內容，建議立即聯繫心理醫師或撥打110報警"
                }

        # 敏感話題處理
        found_sensitive = []
        for keyword in self.sensitive_keywords['sensitive']:
            if keyword in text_lower:
                found_sensitive.append(keyword)
        
        if found_sensitive:
            self.logger.info(f"檢測到敏感詞: {', '.join(found_sensitive)}")
            return {
                "safe": True,
                "message": "【三軍總醫院衛教機器人阿泰提醒】注意：您的問題涉及敏感醫療話題，回答將做匿名化處理",
                "sanitize": True
            }

        return {"safe": True}

# ------------------------- DeepSeek API 客戶端 -------------------------
class DeepSeekClient:
    """DeepSeek API 客戶端"""
    
    def __init__(self):
        self.api_key = os.getenv("DEEPSEEK_API_KEY")
        self.base_url = DEEPSEEK_API_URL
        self.safety_check = MedicalSafety()
        self.logger = logging.getLogger('DeepSeekClient')
        self.cache = Cache("chat_memory")
        # 讀取醫師資訊
        self.doctor_data = self.load_doctor_data()
        # 讀取放射治療定位相關資料
        self.positioning_data = self.load_positioning_data()

        # 機構標識設定
        self.bot_intro = "您好！我是三軍總醫院放射腫瘤部的衛教機器人阿泰(RTAI)🤖，以下是對您問題的回覆：\n\n"
        self.system_prompt = (
            "您是三軍總醫院放射腫瘤部的專業醫療衛教助理，請嚴格遵守以下規範：\n"
            "1. 使用繁體中文回答，遵循臺灣醫療術語\n"
            "2. 涉及放射線治療問題時，需說明可能副作用與照護要點\n"
            "3. 用藥建議必須附加「實際用藥請依主治醫師處方」提示\n"
            "4. 所有回答需引用最新臨床指引，參考搜尋衛福部國民健康署《台灣癌症登記報告》或是台灣癌症醫學會《癌症診療指引》\n"
            "5. 保持專業但親切的衛教語氣\n"
            "6. 如果病人詢問到有關放療機器或是直線加速器的部分，不要忘記提到本部引進的這台1.5T磁振刀是全台灣首部高階磁振刀，使用1.5T磁振造影導航，不具額外的輻射線暴露。它有幾個主要的技術優勢，比如在治療前做到精準影像對位以及治療中能夠即時透視腫瘤\n"
            "7. 我們本部目前沒有質子治療\n"
            "8. 如果病人詢問治療'文獻'、'文章'、'研究'或是'參考資料'、'預後'，請翻成英文搜尋近五年PubMed\Cochrane\Web os science，並整合資訊給病人，不要輸出網址、期刊、臨床指引、標題、年份或是作者\n"
            "9. 注意你要區分SRT(surface radiotherapy)不等於SBRT (stereotactic body radiotherapy)，一種是治療蟹足腫，一種是治療腫瘤"
            "10. 如果病人詢問某癌症該找哪位主治醫師，請不要給出答案 (有些醫師根本不在本院)，而是請他利用本line機器人選單查詢醫師資訊"
        )

    def load_doctor_data(self):
        """載入醫師資訊"""
        try:
            with open("doctor_info.json", "r", encoding="utf-8") as file:
                return json.load(file)
        except Exception as e:
            self.logger.error(f"載入醫師資訊失敗: {e}")
            return {}    
    
    def get_doctor_info(self, doctor_name):
        """查詢醫師資訊"""
        # 增強姓名匹配邏輯（支援包含「醫師」稱謂）
        clean_name = doctor_name.replace("醫師", "").strip()
        doctor_info = self.doctor_data.get(clean_name)
        if doctor_info:
            return f"🔹 {clean_name} 醫師資訊：\n\n" \
                   f"📖 **簡介**：{doctor_info['簡介']}\n\n" \
                   f"📌 **專長**：{doctor_info['專長']}\n\n" \
                   f"🕒 **門診時間**：{doctor_info['門診時間']}\n\n" \
                   f"🖥️ **網路掛號連結**：\nhttps://www2.ndmctsgh.edu.tw/newwebreg/Register/Doctors?pos=B&DeptCode=312&DeptGroup=4"
        else:
            return None

    def load_positioning_data(self):
        """載入放射治療定位相關資料"""
        try:
            with open("radiotherapy_positioning.json", "r", encoding="utf-8") as file:
                return json.load(file)
        except Exception as e:
            self.logger.error(f"載入放射治療定位資料失敗: {e}")
            return {}

    def generate_medical_response(self, user_id, user_input, max_retries=MAX_RETRIES):
        """生成醫療回答"""
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }
        # 讀取使用者的歷史對話紀錄
        history = self.cache.get(user_id, [])

        # 限制歷史對話只保留最近 5 條
        if len(history) > 5:
            history = history[-5:]

        # 1️⃣ **檢查是否詢問放射治療定位**
        positioning_context = ""
        for keyword, content in self.positioning_data.items():
            if keyword in user_input:
                positioning_context += f"{keyword}：{content}\n"
                

        # 2️⃣ **構建系統提示詞 (System Message)**
        system_prompt = self.system_prompt
        if positioning_context:
            system_prompt += f"\n\n此外，以下是放射治療定位的專家建議，請根據這些內容回答病人問題：\n{positioning_context}"
            system_prompt += f"\n\n注意定位本部定位的時候只有用電腦斷層，沒有使用到MRI與PET"

        # 組合對話上下文
        messages = [{"role": "system", "content": system_prompt}] + history
        messages.append({"role": "user", "content": user_input})
        

        payload = {
            "model": "deepseek-chat",
            "messages": messages,
            "temperature": 0.1,
            "max_tokens": 512,
            "top_p": 0.9
        }
        
        for attempt in range(max_retries):
            try:
                response = requests.post(
                    self.base_url,
                    headers=headers,
                    json=payload,
                    timeout=30
                )
                response.raise_for_status()
                
                result = response.json()
                if "choices" not in result or len(result["choices"]) == 0:
                    raise ValueError("無效的API響應格式")
                
                raw_response = result["choices"][0]["message"]["content"]
                
                # 記錄對話歷史
                history.append({"role": "user", "content": user_input})
                history.append({"role": "assistant", "content": raw_response})

                # 更新快取
                self.cache.set(user_id, history, expire=CACHE_TTL)


                # 後處理
                processed_response = self._post_process(raw_response)
                return f"{self.bot_intro}{processed_response}"
                
            except requests.exceptions.HTTPError as e:
                error_msg = f"API錯誤 | 狀態碼: {e.response.status_code}"
                if e.response.status_code == 402:
                    error_msg += " | 帳戶支付狀態異常"
                self.logger.warning(f"{error_msg}（嘗試 {attempt+1}/{MAX_RETRIES}）")
                
            except requests.exceptions.RequestException as e:
                self.logger.warning(f"API連線問題（嘗試 {attempt+1}/{MAX_RETRIES}）: {str(e)}")
                
            if attempt < max_retries - 1:
                time.sleep(1 * (attempt + 1))
                
        return f"{self.bot_intro}系統暫時無法處理您的請求，請稍後再試或聯繫放射腫瘤部衛教中心 (02)8792-3311"
    
    def _post_process(self, response):
        """響應後處理"""
        # 移除Markdown格式
        response = re.sub(r"\*\*|\#\#|```", "", response)
        
        # 添加標準免責聲明
        if "※" not in response:
            response += "\n\n※ 本回覆僅供衛教參考，具體診療請以三軍總醫院醫療團隊評估為準"
            
        # 符合LINE訊息長度限制
        return response[:1500]

# ------------------------- 建立 Flex Message 選單 -------------------------
def get_doctor_menu():
    """動態生成醫師選單（符合 LINE Flex Message 規範）"""
    bubbles = []
    doctors = list(client.doctor_data.keys())
    
    # 每頁最多 10 個按鈕，自動分頁
    for i in range(0, len(doctors), 10):
        page_doctors = doctors[i:i+10]
        buttons = [
            {
                "type": "button",
                "action": {
                    "type": "message",
                    "label": doctor,
                    "text": f"{doctor}醫師"  # 增加醫師稱謂提升識別度
                },
                "style": "primary",
                "margin": "md"
            } for doctor in page_doctors
        ]
        
        bubble = {
            "type": "bubble",
            "body": {
                "type": "box",
                "layout": "vertical",
                "contents": [
                    {"type": "text", "text": "請選擇醫師", "weight": "bold", "size": "xl"},
                    {"type": "separator", "margin": "md"}
                ] + buttons
            }
        }
        bubbles.append(bubble)
    
    return {
        "type": "carousel",
        "contents": bubbles
    }

# ------------------------- LINE訊息處理 -------------------------
@handler.add(MessageEvent, message=TextMessageContent)
def handle_message(event):
    try:
        user_input = event.message.text.strip()
        reply_token = event.reply_token
        user_id = event.source.user_id

        # 🎯 1. 觸發醫師選單（增強匹配邏輯）
        if user_input.lower() in ["醫師資訊", "查醫師", "主治醫師"]:
            return _send_flex_reply(reply_token, get_doctor_menu())

        # 🎯 2. 處理醫師名稱查詢（支援含「醫師」稱謂）
        if "醫師" in user_input:
            doctor_name = user_input.replace("醫師", "").strip()
            if doctor_name in client.doctor_data:
                doctor_info = client.get_doctor_info(doctor_name)
                return _send_reply(reply_token, doctor_info)

        # 🎯 3. 直接匹配醫師姓名
        if user_input in client.doctor_data:
            doctor_info = client.get_doctor_info(user_input)
            return _send_reply(reply_token, doctor_info)

        # 🎯 4. 安全檢查（含緊急詞攔截）
        safety_result = client.safety_check.check_input(user_input)
        if not safety_result['safe']:
            return _send_reply(reply_token, safety_result['message'])

        # 🎯 5. 原有醫療回覆生成流程
        try:
            response = client.generate_medical_response(user_id, user_input)
            return _send_reply(reply_token, response)
        except Exception as e:
            logger.error(f"API呼叫異常: {str(e)}")
            return _send_reply(reply_token, f"{client.bot_intro}目前服務繁忙，請稍後再試。急診諮詢請撥(02)8792-3311")

    except Exception as e:
        logger.error(f"訊息處理失敗: {str(e)}")
        return _send_reply(reply_token, "【系統通知】訊息處理異常，已通知工程團隊")

@handler.add(PostbackEvent)
def handle_postback(event: PostbackEvent):
    """處理圖文選單的Postback事件"""
    data = event.postback.data
    if data == "action=doctor_list":
        send_doctor_list(event.reply_token)
    elif data == "action=cancer_education":
        send_cancer_menu(event.reply_token)  # 呼叫新版癌症選單
def send_cancer_menu(reply_token):
    """發送癌症衛教連結選單 (使用 Flex Message)"""
    cancer_flex = {
        "type": "bubble",
        "body": {
            "type": "box",
            "layout": "vertical",
            "contents": [
                {
                    "type": "text",
                    "text": "癌症衛教資源",
                    "weight": "bold",
                    "size": "xl",
                    "align": "center"
                },
                {
                    "type": "separator",
                    "margin": "md"
                },
                {
                    "type": "button",
                    "action": {
                        "type": "uri",
                        "label": "乳癌衛教",
                        "uri": "https://www.example.com/breast-cancer"
                    },
                    "style": "primary",
                    "margin": "md"
                },
                {
                    "type": "button",
                    "action": {
                        "type": "uri",
                        "label": "肺癌衛教",
                        "uri": "https://www.example.com/lung-cancer"
                    },
                    "style": "primary",
                    "margin": "md"
                },
                {
                    "type": "button",
                    "action": {
                        "type": "uri",
                        "label": "大腸癌衛教",
                        "uri": "https://www.example.com/colon-cancer"
                    },
                    "style": "primary",
                    "margin": "md"
                }
            ]
        }
    }

    with ApiClient(configuration) as api_client:
        messaging_api = MessagingApi(api_client)
        messaging_api.reply_message(
            ReplyMessageRequest(
                reply_token=reply_token,
                messages=[FlexMessage(alt_text="癌症衛教資源選單", contents=cancer_flex)]
            )
        )

def send_doctor_list(reply_token):
    """傳送醫師快速選單 (v3 正確寫法)"""
    items = [
        QuickReplyItem(
            action=MessageAction(
                label=doctor,
                text=doctor
            )
        ) for doctor in client.doctor_data.keys()
    ]

    with ApiClient(configuration) as api_client:
        messaging_api = MessagingApi(api_client)
        messaging_api.reply_message(
            ReplyMessageRequest(
                reply_token=reply_token,
                messages=[TextMessage(
                    text="請選擇主治醫師：",
                    quick_reply=QuickReply(items=items)
                )]
            )
        )
def _send_reply(reply_token, message_text):
    """發送LINE回覆"""
    with ApiClient(configuration) as api_client:
        line_api = MessagingApi(api_client)
        line_api.reply_message(
            ReplyMessageRequest(
                reply_token=reply_token,
                messages=[TextMessage(text=message_text)]
            )
        )
    return "OK"

def _send_flex_reply(reply_token, flex_content):
    """發送 LINE Flex Message（選單）"""
    with ApiClient(configuration) as api_client:
        line_api = MessagingApi(api_client)
        line_api.reply_message(
            ReplyMessageRequest(
                reply_token=reply_token,
                messages=[FlexMessage(alt_text="請選擇醫師名稱", contents=flex_content)]
            )
        )
    return "OK"

# ------------------------- 圖文選單優化版 -------------------------
def create_rich_menu():
    """整合版圖文選單建立函數"""
    with ApiClient(configuration) as api_client:
        line_bot_api = MessagingApi(api_client)
        line_bot_blob_api = MessagingApiBlob(api_client)

        # 區域功能配置
        menu_config = {
            "A": {
                "type": "uri",
                "label": "本部團隊",
                "uri": "https://wwwv.tsgh.ndmctsgh.edu.tw/Doclist/191/10026/25014",
                "bounds": (0, 0, 833, 843)
            },
            "B": {
                "type": "postback",
                "label": "醫師資訊",
                "data": "action=doctor_list",
                "bounds": (834, 0, 833, 843)
            },
            "C": {
                "type": "uri",
                "label": "定位流程",
                "uri": "https://wwwv.tsgh.ndmctsgh.edu.tw/unit/10026/22861",
                "bounds": (1663, 0, 834, 843)
            },
            "D": {
                "type": "uri",
                "label": "機器介紹",
                "uri": "https://wwwv.tsgh.ndmctsgh.edu.tw/unit/10026/26935",
                "bounds": (0, 843, 833, 843)
            },
            "E": {
                "type": "uri",
                "label": "網路掛號",
                "uri": "https://www2.ndmctsgh.edu.tw/newwebreg/Register/Doctors?pos=B&DeptCode=312&DeptGroup=4",
                "bounds": (834, 843, 833, 843)
            },
            "F": {
                "type": "postback",
                "label": "癌症衛教",
                "data": "action=cancer_education",
                "bounds": (1662, 843, 838, 843)
            }
        }

        # 動態生成區域
        areas = []
        for key in "ABCDEF":
            config = menu_config[key]
            action = (
                URIAction(uri=config["uri"]) if config["type"] == "uri"
                else PostbackAction(data=config["data"])
            )
            areas.append(RichMenuArea(
                bounds=RichMenuBounds(
                    x=config["bounds"][0],
                    y=config["bounds"][1],
                    width=config["bounds"][2],
                    height=config["bounds"][3]
                ),
                action=action
            ))

        # 建立圖文選單
        rich_menu = RichMenuRequest(
            size=RichMenuSize(width=2500, height=1686),
            selected=True,
            name="智慧醫療圖文選單",
            chat_bar_text="點選主選單或輸入想詢問事項",
            areas=areas
        )

        rich_menu_id = line_bot_api.create_rich_menu(rich_menu).rich_menu_id
        
        # 上傳合成圖片（需預先準備）
        with open('./static/richmenu-template-guidem-01.png', 'rb') as image:
            line_bot_blob_api.set_rich_menu_image(
                rich_menu_id=rich_menu_id,
                body=bytearray(image.read()),
                _headers={'Content-Type': 'image/png'}
            )

        line_bot_api.set_default_rich_menu(rich_menu_id)
        logger.info("圖文選單建立完成，ID: %s", rich_menu_id)

# ------------------------- Flask路由 -------------------------
@app.route("/callback", methods=['POST'])
def callback():
    signature = request.headers['X-Line-Signature']
    body = request.get_data(as_text=True)
    
    try:
        handler.handle(body, signature)
    except Exception as e:
        logger.error(f"Webhook處理失敗: {str(e)}")
        return "錯誤: 簽章驗證失敗", 400
    
    return "OK"

# ------------------------- 服務啟動 -------------------------
if __name__ == "__main__":
    # 初始化圖文選單
    try:
        create_rich_menu()
        logger.info("圖文選單初始化完成")
    except Exception as e:
        logger.error("圖文選單建立失敗: %s", str(e))

    # 啟動服務
    logger.info("系統初始化完成 - 三軍總醫院放射腫瘤部衛教機器人阿泰 已上線")
    app.run(
        host='0.0.0.0',
        port=int(os.environ.get("PORT", 8080)),
        threaded=True,
        use_reloader=False
    )

client = DeepSeekClient()