from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from Bio import Entrez, Medline
from itertools import count
from typing import List, Optional
import uvicorn
import re
import time
from dotenv import load_dotenv
import os

Entrez.email = os.getenv("ENTREZ_EMAIL")

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],    # dev: allow all; prod: lock this down to your React origin
    allow_methods=["GET"],
    allow_headers=["*"],
)
def strip_tags(text):
  return re.sub(r'<[^>]+>', '', text)

def fetch_pubmed_records( affiliation: str, from_date: str, to_date: Optional[str] ) -> List[dict]:
  if not to_date or to_date.strip() == "":
    to_date = "Present"
  query = f"{affiliation}[AFFL] AND ('{from_date}'[PDAT] : '{to_date}'[PDAT])"

  h = Entrez.esearch(db="pubmed", term=query, retmax=10000, sort="pub+date")
  rec = Entrez.read(h); h.close()

  if not rec["IdList"]:
    return []
  
  id_list = rec["IdList"]
  
  all_articles = []
  for i in range(0, len(id_list), 200):
    ids_chunk = id_list[i:i+200]
    handle = Entrez.efetch(db="pubmed", id=ids_chunk, rettype="xml", retmode="text")
    articles = Entrez.read(handle)
    handle.close()
    all_articles.extend(articles["PubmedArticle"])
    time.sleep(0.5)

  out = []
  ## Iterating over each article and saving in df
  for article in all_articles:
    article_info = article["MedlineCitation"]["Article"]
    pubdata = article["PubmedData"]

    ## Collaborations
    authors_list = article_info.get("AuthorList", [])
    collab = ""
    collab_author = ""
    try:
      collab = strip_tags(article_info.get("AuthorList", [])[-1]["CollectiveName"])
      collab_author = collab
    except:
      if authors_list[-1].get("LastName", "").lower() in ["collaboration", "group", "consortium"]:
        collab = authors_list[-1].get("LastName", "") + ", " + authors_list[-1].get("ForeName", "")
      else:
        collab = ""

    ## Publication Type
    general_types = {"Journal Article"}
    pub_types = [pt for pt in article_info.get("PublicationTypeList", []) if pt not in general_types]
    pub_type = ", ".join(pub_types) if pub_types else "Journal Article"

    ## Article PMID
    pmid = article["MedlineCitation"]["PMID"]

    ## Article Title
    title = strip_tags(article_info["ArticleTitle"])

    ## Journal Title
    journal = strip_tags(article_info["Journal"]["Title"])

    ## Article Language
    lang = strip_tags(article_info["Language"][0])

    ## DOI Number
    doi = ""
    for elink in article_info.get("ELocationID", []):
      if elink.attributes.get("EIdType") == "doi":
        doi = elink

    ## DOI Link
    doi_link = ""
    if doi != "":
      doi_link = f"https://doi.org/{doi}"

    ## Authors - Affiliation
    counter = count(1)
    authors_list = article_info.get("AuthorList", [])
    authors_affiliation = "\n".join([f"{next(counter)}. {a.get('LastName', '')} {a.get('ForeName', '')} - {b.get('Affiliation', '')}" for a in authors_list for b in a["AffiliationInfo"] if "LastName" in a])
    if collab_author != "" and authors_affiliation != "":
      authors_affiliation += f"\n{collab_author}"
    elif collab_author != "":
      authors_affiliation = collab_author

    ## First and Last Author
    first_author = ""
    last_author = ""

    if authors_list:
      try:
        # First Author
        first = authors_list[0]
        first_name = f"{first.get('ForeName', '')} {first.get('LastName', '')}".strip()
        first_affil = first.get("AffiliationInfo", [{}])[0].get("Affiliation", "")
        first_author = ", ".join(filter(None, [first_name, first_affil]))

        # Last Author
        if collab != "":
          last = authors_list[-2]
        else:
          last = authors_list[-1]
        last_name = f"{last.get('ForeName', '')} {last.get('LastName', '')}".strip()
        last_affil = last.get("AffiliationInfo", [{}])[0].get("Affiliation", "")
        last_author = ", ".join(filter(None, [last_name, last_affil]))

      except:
        first_author = ""
        last_author = ""


    ## PUB Link
    pub_link = "https://pubmed.ncbi.nlm.nih.gov/"+str(pmid)

    ## Dates
    pub_date = article_info["Journal"]["JournalIssue"]["PubDate"].get("Year", "N/A")
    month = article_info["Journal"]["JournalIssue"]["PubDate"].get("Month", "")
    article_date = article_info.get("ArticleDate", [])
    if article_date != []:
      article_date = strip_tags("/".join([value for key,value in article_info.get("ArticleDate", [])[0].items()]))

    ## Publication Status
    pub_status = pubdata.get("PublicationStatus", "")

    ## PMCID
    pmcid = ""
    article_ids = pubdata.get("ArticleIdList", [])
    for a in article_ids:
      if a.attributes.get("IdType") == "pmc":
        pmcid = f"PMCID: {a}"

    ## Citation
    # Authors
    author_str = strip_tags(", ".join(
        f"{a.get('LastName', '')} {a.get('Initials', '')}"
        for a in authors_list if "LastName" in a and "Initials" in a
    ))
    if collab_author != "" and author_str != "":
      author_str += f"; {collab_author}"
    elif collab_author != "":
      author_str = collab_author

    # Source
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
    records = Medline.parse(handle)
    record = next(records)  # Only one PMID, so get the first record
    time.sleep(0.1)

    source = record.get("SO", "N/A")

    idx = source.find("eCollection")
    if idx != -1:
      source = source[:idx-1].strip()
    
    if pub_status == "aheadofprint":
      source += " Epub ahead of print."

    # Erratum
    erratum = ""
    for ref in article["MedlineCitation"].get("CommentsCorrectionsList", []):
      if ref.attributes.get("RefType") == "ErratumIn":
        refsource = ref.get("RefSource", "").strip()
        erratum = f"Erratum in: {refsource}"
        break  # Only include the first erratum if there are multiple


    citation = f"{author_str}. {title} {source} {erratum} PMID: {pmid}"
    if pmcid:
      citation += f"; {pmcid}."

    out.append({
      "Pubmed ID":                  pmid,
      "DOI":                        str(doi),
      "Title":                      title,
      "Authors":                    authors_affiliation,
      "Journal":                    journal,
      "Article Date":               article_date,
      "Year":                       pub_date,
      "Month of Publication":       month,
      "Type":                       pub_type,
      "Publication Status":         pub_status,
      "Language":                   lang,
      "Link":                       pub_link,
      "DOI link":                   doi_link,
      "Manuscript":                 citation,
      "First Author":               first_author,
      "Last Author":                last_author,
      "Quartile (Scimago)":         "",
      "Impact Factor":              "",
      "Department":                 "",
      "Sector":                     "",
      "Open Access":                "",
      "KU":                         "",
      "MOH":                        "",
      "Other Nation Institutions":  "",
      "National All":               "",
      "Regional":                   "",
      "International":              "",
      "Total International":        ""
    })
  return out

@app.get("/articles")
def articles_endpoint(affiliation: str, from_date: str, to_date: Optional[str] = ""):
  handle = Entrez.einfo(db="pubmed")
  pubmed_db = Entrez.read(handle)
  handle.close()

  db_name = pubmed_db["DbInfo"]["DbName"]
  description = pubmed_db["DbInfo"]["Description"]
  count = pubmed_db["DbInfo"]["Count"]
  last_update = pubmed_db["DbInfo"]["LastUpdate"]

  db_info = {
    "Name": db_name,
    "Description": description,
    "Count": count,
    "Last Update": last_update
  }

  results = fetch_pubmed_records(affiliation, from_date, to_date)
  return {"count": len(results), "articles": results, "db_info": db_info}