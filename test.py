from Bio import Entrez
import json
import time
from googletrans import Translator  # 新增此行
import re
import logging
logger = logging.getLogger(__name__)

def translate_to_english(text):
    """使用 Google Translate 將中文轉換為英文"""
    translator = Translator()
    translation = translator.translate(text, src='zh-CN', dest='en')
    return translation.text

def search_pubmed(keyword, max_results=3):
    """PubMed文獻搜索，支援中文關鍵字轉英文"""
    Entrez.email = "he165076373@hotmail.com"  # 需申請NCBI帳號

    # 🔹 如果輸入為中文，先翻譯成英文
    if re.search("[\u4e00-\u9fff]", keyword):  # 檢測是否包含中文
        keyword = translate_to_english(keyword)
        logger.info(f"🔄 已將關鍵字翻譯為英文: {keyword}")

    handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_results, sort="relevance")
    result = Entrez.read(handle)
    handle.close()
    return result.get("IdList", [])

def fetch_article_details(pubmed_id):
    """安全解析 PubMed 文獻詳情"""
    try:
        # 🔹 明確指定返回格式為 XML
        handle = Entrez.efetch(
            db="pubmed",
            id=pubmed_id,
            retmode="xml",  # 必須指定 XML 格式
            rettype="medline"
        )
        
        # 🔹 分層解析，避免直接訪問鍵值
        articles = Entrez.read(handle)
        if not articles:
            print(f"警告：PMID {pubmed_id} 無有效數據")
            return None
        
        medline_citation = articles[0].get("MedlineCitation", {})
        article = medline_citation.get("Article", {})
        
        # 🔹 處理標題
        title = article.get("ArticleTitle", "無標題資訊")
        
        # 🔹 處理作者列表
        authors = []
        author_list = article.get("AuthorList", [])
        for author in author_list:
            # 處理姓名缺失情況
            last_name = author.get("LastName", "Unknown")
            initials = author.get("Initials", "")
            authors.append(f"{last_name} {initials}")
        
        # 🔹 返回結構化數據
        return {
            "title": title,
            "authors": authors[:3],  # 最多取前3位作者
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
        }
        
    except Exception as e:
        # 🔹 打印具體錯誤訊息
        print(f"解析 PMID {pubmed_id} 失敗，錯誤原因: {str(e)}")
        return None

def handle_message(user_input):
    
        # 🎯 3. 文獻請求觸發機制（新增部分）
        literature_triggers = ["文獻", "研究", "來源", "source", "reference"]
        if any(trigger in user_input.lower() for trigger in literature_triggers):
            print('If any')
            try:
                # PubMed API調用
                pubmed_ids = search_pubmed(user_input, max_results=5)
                print(pubmed_ids)
                articles = [fetch_article_details(pid) for pid in pubmed_ids]
                print(articles)
                
                if not articles:
                    response = "⚠️ 目前未找到相關文獻，建議簡化關鍵字或諮詢醫師。"
                else:
                    response = "📚 以下為PubMed文獻：\n\n"
                    for art in articles:
                        response += f"► {art['title']}\n作者：{', '.join(art['authors'][:2])}\n連結：{art['url']}\n\n"
                    response += "※ 注意：此為學術資料，具體診療請遵醫囑"
                    
                return response
            
            except Exception as e:
                return "文獻服務暫時不可用，請稍後再試"

a = handle_message('食道癌放療劑量提升文獻')
print(a)