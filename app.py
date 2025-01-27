有辦法醫師諮詢

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
from linebot.v3.messaging import Configuration, ApiClient, MessagingApi, ReplyMessageRequest, TextMessage, FlexMessage
from linebot.v3.webhooks import MessageEvent, TextMessageContent
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

        # 🧬 新增醫療服務配置
        self.registration_keywords = ["挂号", "掛號", "預約", "門診", "看診"]
        self.title_suffixes = ["医师", "医生", "醫師", "大夫", "主任", "醫生"]
        self.alias_to_doctor = {}  
        for doctor, info in self.doctor_data.items():
            for alias in info["aliases"]:
                self.alias_to_doctor[alias.lower()] = doctor




        # 讀取資料
        self.doctor_data = self.load_doctor_data()  # 🧬 修改後的載入方法
        self.positioning_data = self.load_positioning_data()




        # 機構標識設定
        self.bot_intro = "您好！我是三軍總醫院放射腫瘤部的衛教機器人阿泰(RTAI)🤖，以下是對您問題的回覆：\n\n"
        self.system_prompt = self.system_prompt = (
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
            )

    # 🧬 修改後的醫師資料載入方法
    def load_doctor_data(self):
        """載入醫師資訊並生成別名"""
        try:
            with open("doctor_info.json", "r", encoding="utf-8") as file:
                doctor_data = json.load(file)
                
                for name, info in doctor_data.items():
                    full_name = str(name).strip()
                    aliases = []
                    
                    # 生成識別別名
                    parts = full_name.split()
                    if parts:
                        surname = parts[0][0]  # 取姓氏第一個字
                        for suffix in self.title_suffixes:
                            aliases.append(f"{surname}{suffix}")
                            aliases.append(f"{full_name}{suffix}")
                            
                    aliases.append(full_name)
                    info["aliases"] = list(set(aliases))
                
                return doctor_data
        except Exception as e:
            self.logger.error(f"載入醫師資訊失敗: {e}")
            return {}

    def load_positioning_data(self):
        """載入放射治療定位相關資料"""
        try:
            with open("radiotherapy_positioning.json", "r", encoding="utf-8") as file:
                return json.load(file)
        except Exception as e:
            self.logger.error(f"載入放射治療定位資料失敗: {e}")
            return {}

    # 🧬 新增醫療特徵檢測方法
    def detect_medical_mentions(self, user_input):
        input_lower = user_input.lower()
        found_doctors = set()
    # 使用集合操作优化匹配
        for alias, doctor in self.alias_to_doctor.items():
            if alias in input_lower:
                found_doctors.add(doctor)
        return {"doctors": list(found_doctors)}

    # 🧬 新增醫療資訊構建方法
    def build_medical_response(self, detection_result):
        """構建醫療相關回應"""
        response = ""
        
        if detection_result["doctors"]:
            response += "🏥 **相關醫師資訊**：\n\n"
            for doctor in list(set(detection_result["doctors"])):  # 去重
                info = self.get_doctor_info(doctor)
                if info:
                    response += info + "\n\n"
                    
        if detection_result["needs_registration"]:
            response += "📅 **掛號服務**：\n"
            response += "1. 放射腫瘤部網路掛號系統：\n"
            response += "   https://www2.ndmctsgh.edu.tw/newwebreg/Register/Doctors?pos=B&DeptCode=312&DeptGroup=4\n"
            response += "2. 人工掛號專線：(02)8792-3311 轉分機 12345\n"
            response += "3. 現場掛號：門診大樓1號櫃台\n\n"
            
        return response

    def get_doctor_info(self, doctor_name):
        """查詢醫師資訊"""
        doctor_info = self.doctor_data.get(doctor_name)
        if doctor_info:
            return f"🔹 {doctor_name} 醫師資訊：\n\n" \
                   f"📖 **簡介**：{doctor_info['簡介']}\n\n" \
                   f"📌 **專長**：{doctor_info['專長']}\n\n" \
                   f"🕒 **門診時間**：{doctor_info['門診時間']}\n\n" \
                   f"🖥️ **網路掛號連結**：\nhttps://www2.ndmctsgh.edu.tw/newwebreg/Register/Doctors?pos=B&DeptCode=312&DeptGroup=4"
        else:
            return None

    def generate_medical_response(self, user_id, user_input, max_retries=MAX_RETRIES):
        """生成醫療回答"""
        # 🧬 新增前置檢測流程
        # 1. 安全檢查
        if not self.safety_check.validate(user_input):
            return self.bot_intro + "您的問題涉及專業醫療建議，建議直接諮詢主治醫師。"
            
        # 2. 快取檢查
        cached = self.cache.get(user_input)
        if cached:
            return self.bot_intro + cached
            
        # 3. 醫療特徵檢測
        detection = self.detect_medical_mentions(user_input)
        medical_response = self.build_medical_response(detection)
        if medical_response:
            final_response = self.bot_intro + medical_response
            self.cache.set(user_input, final_response)
            return final_response

        # 以下為原有API調用流程
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }
        history = self.cache.get(user_id, [])
        if len(history) > 5:
            history = history[-5:]

        positioning_context = ""
        for keyword, content in self.positioning_data.items():
            if keyword in user_input:
                positioning_context += f"{keyword}：{content}\n"

        system_prompt = self.system_prompt
        if positioning_context:
            system_prompt += f"\n\n此外，以下是放射治療定位的專家建議，請根據這些內容回答病人問題：\n{positioning_context}"
            system_prompt += f"\n\n注意定位本部定位的時候只有用電腦斷層，沒有使用到MRI與PET"

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
                
                history.append({"role": "user", "content": user_input})
                history.append({"role": "assistant", "content": raw_response})
                self.cache.set(user_id, history, expire=CACHE_TTL)

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
        response = re.sub(r"\*\*|\#\#|```", "", response)
        if "※" not in response:
            response += "\n\n※ 本回覆僅供衛教參考，具體診療請以三軍總醫院醫療團隊評估為準"
        return response[:1500]

# ------------------------- 建立 Flex Message 選單 -------------------------
def get_doctor_menu():
    """動態生成醫師選單（符合 LINE Flex Message 規範）"""
    bubbles = []
    doctors = list(client.doctor_data.keys())
    
    for i in range(0, len(doctors), 10):
        page_doctors = doctors[i:i+10]
        buttons = [
            {
                "type": "button",
                "action": {
                    "type": "message",
                    "label": doctor,
                    "text": doctor
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

        # 🧬 增強醫療特徵檢測
        detection = client.detect_medical_mentions(user_input)
        if detection["doctors"] or detection["needs_registration"]:
            response = client.build_medical_response(detection)
            return _send_reply(reply_token, client.bot_intro + response)

        # 原有處理流程
        if user_input == "我想查詢我的放射治療主治醫師":
            return _send_flex_reply(reply_token, get_doctor_menu())
            
        if user_input in client.doctor_data:
            doctor_info = client.get_doctor_info(user_input)
            return _send_reply(reply_token, doctor_info)

        safety_result = client.safety_check.check_input(user_input)
        if not safety_result['safe']:
            return _send_reply(reply_token, safety_result['message'])

        try:
            response = client.generate_medical_response(user_id, user_input)
            return _send_reply(reply_token, response)
        except Exception as e:
            logger.error(f"API呼叫異常: {str(e)}")
            return _send_reply(reply_token, f"{client.bot_intro}目前服務繁忙，請稍後再試。急診諮詢請撥(02)8792-3311")

    except Exception as e:
        logger.error(f"訊息處理失敗: {str(e)}")
        return _send_reply(reply_token, "【系統通知】訊息處理異常，已通知工程團隊")

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

client = DeepSeekClient()
# ------------------------- 服務啟動 -------------------------
if __name__ == "__main__":
    logger.info("系統初始化完成 - 三軍總醫院放射腫瘤部衛教機器人阿泰 已上線")
    app.run(
        host='0.0.0.0',
        port=int(os.environ.get("PORT", 8080)),
        threaded=True,
        use_reloader=False
    )
