from wordcloud import *
from database import session, RSIDCitation, Article, article_noun_mapping
from settings import WORDCLOUD_STORAGE
import matplotlib.pyplot as plt
import hashlib
import os


def word_blob(rsid_list):
    noun_list = []
    articles = session.query(Article).join(
            RSIDCitation,
            Article.pmid == RSIDCitation.pmid
        ).filter(
            RSIDCitation.rsid.in_(rsid_list)
        ).distinct().all()
    for pmid in [str(article.pmid) for article in articles]:
        noun_list += article_noun_mapping[pmid]
    return " ".join(noun_list)


# Normalize word counts by number of articles associated with a given RSID, scale by 100
def normalized_word_blob(rsid_list):
    word_counts = {}
    results = session.query(Article, RSIDCitation).join(
            RSIDCitation,
            Article.pmid == RSIDCitation.pmid
        ).filter(
            RSIDCitation.rsid.in_(rsid_list)
        ).all()
    for rsid in rsid_list:
        pmids = [str(result.Article.pmid) for result in results if result.RSIDCitation.rsid == int(rsid)]
        if pmids:
            normalization_factor = 1000/len(pmids)
            rsid_nouns = []
            for pmid in pmids:
                rsid_nouns += article_noun_mapping[pmid]
            for noun in set(rsid_nouns):
                normalized_noun_count = normalization_factor * rsid_nouns.count(noun)
                if noun in word_counts:
                    word_counts[noun] += normalized_noun_count
                else:
                    word_counts[noun] = normalized_noun_count
    normalized_blob = ""
    for word in word_counts:
        normalized_blob += "{} ".format(word) * int(word_counts[word])
    return normalized_blob


# Generate a wordcloud png file from a list of rsids
def create_wordcloud_from_rsids(rsid_list, output_path=None):
    #full_text = word_blob(rsid_list)
    full_text = normalized_word_blob(rsid_list)
    word_cloud = WordCloud(width=1600, height=800).generate(full_text)
    plt.figure(figsize=(20, 10), facecolor='k')
    plt.imshow(word_cloud)
    plt.axis("off")
    plt.tight_layout(pad=0)
    plt.show()
    plt.savefig(output_path, facecolor='k', bbox_inches='tight')


# Function results in IndexError if no results are found
def generate_wordcloud(rsid_list):
    rsid_string = " ".join(rsid_list)
    hash_object = hashlib.sha1(rsid_string.encode())
    hashed_file_name = "{}.png".format(hash_object.hexdigest())
    hashed_file_path = os.path.join(WORDCLOUD_STORAGE, hashed_file_name)
    try:
        create_wordcloud_from_rsids(rsid_list, hashed_file_path)
        return hashed_file_name
    except IndexError:
        return ""
"""
def generate_wordcloud(rsid_list):
    rsid_string = " ".join(rsid_list)
    fleetingcloud = tempfile.SpooledTemporaryFile(max_size=0, mode='w+b', buffering=None, encoding=None, newline=None, suffix='', prefix='tmp', dir=None)
    try:
        create_wordcloud_from_rsids(rsid_list, fleetingcloud)
        return fleetingcloud
    except IndexError:
        return ""
"""
