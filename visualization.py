from wordcloud import *
import matplotlib.pyplot as plt
from database import session, RSIDCitation, Article


def word_blob(rsid_list):
    noun_list = []
    pmid_tuples = session.query(RSIDCitation.pmid).filter(RSIDCitation.rsid.in_(rsid_list)).all()
    pmid_list = [pmid_tuple[0] for pmid_tuple in pmid_tuples]
    articles = session.query(Article).filter(Article.pmid.in_(pmid_list)).all()
    for article in articles:
        noun_list += article.abstract_nouns()
    return " ".join(noun_list)


# Generate a wordcloud png file from a list of rsids
def wordcloud_from_rsids(rsid_list, output_path=None):
    noun_list = []
    pmid_tuples = session.query(RSIDCitation.pmid).filter(RSIDCitation.rsid.in_(rsid_list)).all()
    pmid_list = [pmid_tuple[0] for pmid_tuple in pmid_tuples]
    articles = session.query(Article).filter(Article.pmid.in_(pmid_list)).all()
    for article in articles:
        noun_list += article.abstract_nouns()
    full_text = " ".join(noun_list)
    word_cloud = WordCloud().generate(full_text)
    if output_path is not None:
        plt.imshow(word_cloud)
        plt.axis("off")
        plt.savefig(output_path)
    return word_cloud
