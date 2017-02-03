from wordcloud import *
from database import session, RSIDCitation, Article
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
    for article in articles:
        noun_list += article.abstract_nouns()
    return " ".join(noun_list)


# Generate a wordcloud png file from a list of rsids
def create_wordcloud_from_rsids(rsid_list, output_path=None):
    full_text = word_blob(rsid_list)
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
