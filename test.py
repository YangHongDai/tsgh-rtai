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
    """獲取文獻詳情（需自訂解析邏輯）"""
    try:
        handle = Entrez.efetch(db="pubmed", id=pubmed_id, retmode="xml")
        article_data = Entrez.read(handle)[0]['MedlineCitation']
        
        authors = [f"{author['LastName']} {author['Initials']}" 
                  for author in article_data.get('Article', {}).get('AuthorList', [])]
        
        return {
            'title': article_data['Article']['ArticleTitle'],
            'authors': authors[:3],  # 最多取3位作者
            'url': f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
        }
    except Exception as e:
        logger.error(f"文獻解析失敗 {pubmed_id}: {str(e)}")
        return None

a = search_pubmed("前列腺癌 放射治療", max_results=5)
print(a)