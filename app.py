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

# ------------------------- DeepSeek API 客戶端 -------------------------
class DeepSeekClient:
    """DeepSeek API 客戶端"""
    
    def __init__(self):
        self.api_key = os.getenv("DEEPSEEK_API_KEY")
        self.base_url = DEEPSEEK_API_URL
        self.cache = Cache("chat_memory")
        self.logger = logging.getLogger('DeepSeekClient')
        self.doctor_data = self.load_doctor_data()
        self.positioning_data = self.load_positioning_data()

    def load_doctor_data(self):
        """載入醫師資訊"""
        try:
            with open("doctor_info.json", "r", encoding="utf-8") as file:
                return json.load(file)
        except Exception as e:
            self.logger.error(f"載入醫師資訊失敗: {e}")
            return {}    

    def load_positioning_data(self):
        """載入放射治療定位資料"""
        try:
            with open("radiotherapy_positioning.json", "r", encoding="utf-8") as file:
                return json.load(file)
        except Exception as e:
            self.logger.error(f"載入放射治療定位資料失敗: {e}")
            return {}

    def get_doctor_info(self, doctor_name):
        """查詢醫師資訊"""
        doctor_info = self.doctor_data.get(doctor_name)
        if doctor_info:
            return f"🔹 {doctor_name} 醫師資訊：\n\n" \
                   f"📖 **簡介**：{doctor_info['簡介']}\n\n" \
                   f"📌 **專長**：{doctor_info['專長']}\n\n" \
                   f"🕒 **門診時間**：{doctor_info['門診時間']}\n\n" \
                   f"🖥️ **網路掛號連結**：\nhttps://www2.ndmctsgh.edu.tw/newwebreg/Register"
        else:
            return "查無此醫師資訊，請確認姓名是否正確。"

    def generate_medical_response(self, user_id, user_input, max_retries=MAX_RETRIES):
        """生成醫療回答"""
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }

        history = self.cache.get(user_id, [])
        if len(history) > 5:
            history = history[-5:]

        messages = [{"role": "system", "content": "請根據醫療專業回答使用者的問題"}] + history
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
                response = requests.post(self.base_url, headers=headers, json=payload, timeout=30)
                response.raise_for_status()
                result = response.json()

                if "choices" not in result or len(result["choices"]) == 0:
                    raise ValueError("無效的API響應格式")

                raw_response = result["choices"][0]["message"]["content"]

                history.append({"role": "user", "content": user_input})
                history.append({"role": "assistant", "content": raw_response})
                self.cache.set(user_id, history, expire=CACHE_TTL)

                return raw_response
            except requests.exceptions.RequestException as e:
                self.logger.warning(f"API連線問題（嘗試 {attempt+1}/{MAX_RETRIES}）: {str(e)}")

        return "系統暫時無法處理您的請求，請稍後再試。"

# ✅ **初始化 DeepSeek API 客戶端**
client = DeepSeekClient()

# ------------------------- 建立 Flex Message 選單 -------------------------
def get_flex_menu():
    return {
        "type": "flex",
        "altText": "請選擇您要諮詢的項目",
        "contents": {
            "type": "bubble",
            "body": {
                "type": "box",
                "layout": "vertical",
                "contents": [
                    {"type": "text", "text": "請選擇您要諮詢的類別", "weight": "bold", "size": "lg"},
                    {"type": "separator"},
                    {
                        "type": "button",
                        "action": {"type": "message", "label": "放射治療副作用", "text": "放射治療副作用"},
                        "style": "primary"
                    },
                    {
                        "type": "button",
                        "action": {"type": "message", "label": "放療技術與設備", "text": "放療技術與設備"},
                        "style": "primary"
                    },
                    {
                        "type": "button",
                        "action": {"type": "message", "label": "預約與門診", "text": "預約與門診"},
                        "style": "primary"
                    }
                ]
            }
        }
    }

# ------------------------- LINE Bot 訊息處理 -------------------------
@handler.add(MessageEvent, message=TextMessageContent)
def handle_message(event):
    try:
        user_input = event.message.text.strip()
        reply_token = event.reply_token
        user_id = event.source.user_id

        if user_input in ["我要諮詢", "諮詢"]:
            return _send_reply(reply_token, get_flex_menu())

        if user_input in client.doctor_data:
            return _send_reply(reply_token, client.get_doctor_info(user_input))

        predefined_responses = {
            "放射治療副作用": "可能會有疲倦、皮膚變紅等副作用，請諮詢主治醫師。",
            "放療技術與設備": "我們使用 1.5T 磁振刀，提供精準影像導航。",
            "預約與門診": "請透過三軍總醫院掛號系統預約門診：https://www2.ndmctsgh.edu.tw/newwebreg/Register"
        }

        response = predefined_responses.get(user_input, client.generate_medical_response(user_id, user_input))

        return _send_reply(reply_token, response)

    except Exception as e:
        logger.error(f"訊息處理失敗: {str(e)}")
        return _send_reply(reply_token, "【系統通知】訊息處理異常，已通知工程團隊")

# ------------------------- 發送 LINE 訊息 -------------------------
def _send_reply(reply_token, message_text):
    with ApiClient(configuration) as api_client:
        line_api = MessagingApi(api_client)
        line_api.reply_message(
            ReplyMessageRequest(
                reply_token=reply_token,
                messages=[TextMessage(text=message_text)]
            )
        )
    return "OK"

# ------------------------- Flask 伺服器啟動 -------------------------
if __name__ == "__main__":
    app.run(host='0.0.0.0', port=int(os.environ.get("PORT", 8080)), threaded=True, use_reloader=False)
