from Bio import Entrez
import xml.etree.ElementTree as ET
from settings import EMAIL
from datetime import date, datetime
from http.client import IncompleteRead


Entrez.email = EMAIL


# Author is initialized with the Author xml element
class Author:
    def __init__(self, data):
        # Collective vs. Person
        if data.find('CollectiveName') is not None:
            self.collective_name = data.find('CollectiveName').text
        else:
            try:
                self.last_name = data.find('LastName').text
            except AttributeError:
                self.last_name = ''
            try:
                self.first_name = data.find('ForeName').text
            except AttributeError:
                self.first_name = ''
            try:
                self.initials = data.find('Initials').text
            except AttributeError:
                self.initials = ''
        self.affiliations = []
        for affiliation in data.iter('Affiliation'):
            self.affiliations.append(affiliation.text)


# Data passed should be PubmedArticle xml tree
class PubmedArticle:
    def __init__(self, data):
        self.abstract = ""
        medline_element = data.find('MedlineCitation')
        self.pmid = int(medline_element.find('PMID').text)
        self.abstract = " ".join([abstract.text for abstract in data.iter('AbstractText') if abstract.text is not None])
        self.authors = [Author(author_element) for author_element in data.iter('Author')]
        date_created_element = medline_element.find('DateCreated')
        date_revised_element = medline_element.find('DateRevised')
        self.date_created = date(
            year=int(date_created_element.find('Year').text),
            month=int(date_created_element.find('Month').text),
            day=int(date_created_element.find('Day').text)
        )
        self.date_revised = date(
            year=int(date_revised_element.find('Year').text),
            month=int(date_revised_element.find('Month').text),
            day=int(date_revised_element.find('Day').text)
        )

    # This method currently has no error-handling, may need to write it to handle timeouts, etc.
    def rsids(self):
        data = Entrez.read(
            Entrez.elink(
                dbfrom='pubmed',
                db='snp',
                linkname='pubmed_snp_cited',
                id=self.pmid
            )
        )
        rsid_list = [id_dict['Id'] for id_dict in data[0]['LinkSetDb'][0]['Link']]
        return rsid_list


# RSID class contains methods to gather associated PubmedArticle objects
class RSID:
    def __init__(self, rsid_number, annotate=False):
        self.rsid_number = rsid_number
        if annotate:
            # Can add additional attributes to mapping below
            attribute_mapping = {
                "CLINICAL_SIGNIFICANCE": "clinical_significance",
                "FXN_CLASS": "fxn_class",
                "VALIDATED": "validated"
            }
            summary_handle = Entrez.esummary(
                db="snp",
                id=self.rsid_number
            )
            summary_tree = ET.fromstring(summary_handle.read())
            summary_handle.close()
            for item in summary_tree.iter('Item'):
                item_attribute = item.attrib
                name_attribute = item_attribute["Name"]
                if item_attribute["Name"] in attribute_mapping:
                    setattr(self, attribute_mapping[name_attribute], item.text)

    def pubmed_articles(self):
        articles = []
        data = Entrez.read(
            Entrez.elink(
                dbfrom='snp',
                db='pubmed',
                linkname='snp_pubmed_cited',
                id=self.rsid_number
            )
        )
        pmid_list = [id_dict['Id'] for id_dict in data[0]['LinkSetDb'][0]["Link"]]
        # Fetch pubmed xml in batches of 10,000
        start = 0
        while start < len(pmid_list):
            pmid_string = ",".join(pmid_list[start:start+10000])
            # Possible to do 4-5 requests concurrently?
            pubmed_fetch_handle = Entrez.efetch(
                db="pubmed",
                rettype="xml",
                retmode="xml",
                retstart=start,
                retmax=10000,
                id=pmid_string,
            )
            pubmed_set_xml = ET.fromstring(pubmed_fetch_handle.read())
            pubmed_fetch_handle.close()
            # should try to find a way to thread this part, might be the bottleneck
            articles += [PubmedArticle(article_xml) for article_xml in pubmed_set_xml.findall('PubmedArticle')]
            start += 10000
        return articles


# Get all rsid objects from snp_pubmed_cited search
def get_all_rsids():
    rsid_list = []
    retstart = 0
    search_string = "snp_pubmed_cited[sb]"
    search_results = Entrez.read(
        Entrez.esearch(
            term=search_string,
            db='snp',
            retstart=retstart,
            retmax=100000,
            usehistory='y'
        )
    )
    number_of_results = search_results["Count"]
    rsid_list += search_results["IdList"]
    additional_queries = int(number_of_results) // 100000
    while additional_queries != 0:
        while True:
            try:
                retstart += 100000
                search_results = Entrez.read(
                    Entrez.esearch(
                        term=search_string,
                        db='snp',
                        retstart=retstart,
                        retmax=100000,
                        usehistory='y'
                    )
                )
                rsid_list += search_results["IdList"]
                additional_queries -= 1
            except IncompleteRead:
                continue
            except TimeoutError:
                continue
            break
    return [RSID(rsid_number) for rsid_number in rsid_list]


# Get pubmed articles, optionally takes in a datetime.date object and returns articles since then, otherwise gets all
def get_pubmed_articles(since_date=None):
    articles = []
    if since_date is not None:
        today = date.today()
        number_of_days = (today - since_date).days
        pubmed_search = Entrez.read(
            Entrez.esearch(
                db='pubmed',
                term='pubmed_snp_cited[sb]',
                retstart=0,
                retmax=100000,
                usehistory='y',
                datetype='mdat',
                reldate=number_of_days,
            )
        )
    else:
        pubmed_search = Entrez.read(
            Entrez.esearch(
                db='pubmed',
                term='pubmed_snp_cited[sb]',
                retstart=0,
                retmax=100000,
                usehistory='y'
            )
        )
    start = 0
    print(pubmed_search["Count"])
    query_key = pubmed_search["QueryKey"]
    web_env=pubmed_search["WebEnv"]
    while start < int(pubmed_search["Count"]):
        while True:
            try:
                print("fetching article results from {} up".format(start))
                n1 = datetime.now()
                pubmed_fetch_handle = Entrez.efetch(
                    db="pubmed",
                    rettype="xml",
                    retmode="xml",
                    retstart=start,
                    retmax=10000,
                    WebEnv=web_env,
                    query_key=query_key,
                    usehistory="y"
                )
                pubmed_set_xml = ET.fromstring(pubmed_fetch_handle.read())
                pubmed_fetch_handle.close()
                articles += [PubmedArticle(article_xml) for article_xml in pubmed_set_xml.findall('PubmedArticle')]
                n2 = datetime.now()
                print("Got articles in {} seconds".format((n2-n1).total_seconds()))
                start += 10000
            except IncompleteRead:
                continue
            except TimeoutError:
                continue
            break
    return articles


# Returns article object(s) from pmid string (one pmid, or comma-seperated list)
def get_pubmed_article_from_pmid(pmid):
    pubmed_fetch_handle = Entrez.efetch(
        db="pubmed",
        rettype="xml",
        retmode="xml",
        retstart=0,
        retmax=10000,
        id=str(pmid),
    )
    pubmed_xml = ET.fromstring(pubmed_fetch_handle.read())
    pubmed_fetch_handle.close()
    return [PubmedArticle(article_xml) for article_xml in pubmed_xml.findall('PubmedArticle')]
