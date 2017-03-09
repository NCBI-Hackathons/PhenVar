from wordcloud import *
from database import session, RSIDCitation, Article, article_noun_mapping
from settings import WORDCLOUD_STORAGE
import matplotlib.pyplot as plt
import hashlib
import os


# Create table - NOUN|RSID's|#Articles
def word_statistics(rsid_list, weights):#, normalization_type):
    # dictionary in form {"noun": [word_count, article_count, rsids]}
    word_stats = {}
    #noun_list = []
    results = session.query(Article, RSIDCitation).join(
        RSIDCitation,
        Article.pmid == RSIDCitation.pmid
    ).filter(
        RSIDCitation.rsid.in_(rsid_list)
    ).all()
    recorded_pmids = []
    for result in results:
        pmid = str(result.Article.pmid)
        rsid = result.RSIDCitation.rsid
        nouns = article_noun_mapping[pmid]
        counting = pmid not in recorded_pmids
        for noun in set(nouns):
            if counting:
                if noun not in word_stats:
                    word_stats[noun] = [0, 0, []]
                word_stats[noun][0] += nouns.count(noun) # Add to word count
                word_stats[noun][1] += 1                 # Add to article count
            if rsid not in word_stats[noun][2]:
                word_stats[noun][2].append(rsid)
        if counting:
            recorded_pmids.append(pmid)
    for noun in word_stats:
        weight = 1
        noun_stats = word_stats[noun]
        weights_dictionary = {
            "word_count": noun_stats[0],
            "article_count": noun_stats[1],
            "rsid_count": len(noun_stats[2]),
        }
        for weight_label in weights:
            weight *= weights_dictionary[weight_label]
        word_stats[noun].append(weight)
    return word_stats


def word_blob(stats):
    blob = ""
    for noun in stats:
        blob += "{} ".format(noun) * stats[noun][3]
    return blob

#def weighted_word_blob(rsid_list, weights):
#    stats = word_statistics(rsid_list)
#    blob = ""
#    for noun in stats:
#        weight = 1
#        noun_stats = stats[noun]
#        weights_dictionary = {
#            "word_count": noun_stats[0],
#            "article_count": noun_stats[1],
#            "rsid_count": len(noun_stats[2]),
#        }
#        for weight_label in weights:
#            weight *= weights_dictionary[weight_label]
#        blob += "{} ".format(noun) * weight
#    return blob



# Creates a blob of words from all articles associated with rsids in a list
#def word_blob(rsid_list):
#    #noun_list = []
#    #articles = session.query(Article).join(
#    #        RSIDCitation,
#    #        Article.pmid == RSIDCitation.pmid
#    #    ).filter(
#    #        RSIDCitation.rsid.in_(rsid_list)
#    #    ).distinct().all()
#    #for pmid in [str(article.pmid) for article in articles]:
#    #    noun_list += article_noun_mapping[pmid]
#    #return " ".join(noun_list)
#    stats = word_statistics(rsid_list)
#    #return " ".join((" ".format(noun) * stats[noun][0] for noun in stats))
#    blob = ""
#    for noun in stats:
#        blob += "{} ".format(noun) * stats[noun][0]
#    return blob


# Normalize word counts by number of articles associated with a given RSID, scale by 1000
def normalized_word_blob_by_rsid(rsid_list):
    stats = word_statistics(rsid_list)
    blob = ""
    for rsid in rsid_list:
        number_of_citations = session.query(RSIDCitation).filter_by(rsid=int(rsid)).count()

    #word_counts = {}
    #results = session.query(Article, RSIDCitation).join(
    #        RSIDCitation,
    #        Article.pmid == RSIDCitation.pmid
    #    ).filter(
    #        RSIDCitation.rsid.in_(rsid_list)
    #    ).all()
    #for rsid in rsid_list:
    #    pmids = [str(result.Article.pmid) for result in results if result.RSIDCitation.rsid == int(rsid)]
    #    if pmids:
    #        normalization_factor = 1000/len(pmids)
    #        rsid_nouns = []
    #        for pmid in pmids:
    #            rsid_nouns += article_noun_mapping[pmid]
    #        for noun in set(rsid_nouns):
    #            normalized_noun_count = normalization_factor * rsid_nouns.count(noun)
    #            if noun in word_counts:
    #                word_counts[noun] += normalized_noun_count
    #            else:
    #                word_counts[noun] = normalized_noun_count
    #normalized_blob = ""
    #for word in word_counts:
    #    normalized_blob += "{} ".format(word) * int(word_counts[word])
    #return normalized_blob


# Normalize word counts by number of articles each word occurs in
def normalized_word_blob_by_articles(rsid_list):
    stats = word_statistics(rsid_list)
    #return " ".join(("{} ".format(noun) * stats[noun][0] * stats[noun][1] for noun in stats))
    blob = ""
    for noun in stats:
        blob += "{} ".format(noun) * stats[noun][0] * stats[noun][1] # noun * word_count * article_count
    return blob
    #word_counts = {}
    #noun_dict = {}
    #noun_list = []
    #pmid_list = (
    #    str(article.pmid) for article in session.query(Article).join(
    #        RSIDCitation,
    #        Article.pmid == RSIDCitation.pmid
    #    ).filter(
    #        RSIDCitation.rsid.in_(rsid_list)
    #    ).distinct().all()
    #)
    #for pmid in pmid_list:
    #    pmid_noun_list = article_noun_mapping[pmid]
    #    noun_dict[pmid] = pmid_noun_list
    #    noun_list += pmid_noun_list
    #for noun in set(noun_list):
    #    total_count = noun_list.count(noun)
    #    number_of_articles = len([pmid for pmid in noun_dict if noun in noun_dict[pmid]])
    #    word_counts[noun] = total_count * number_of_articles
    #normalized_blob = ""
    #for word in word_counts:
    #    normalized_blob += "{} ".format(word) * word_counts[word]
    #return normalized_blob


# Create wordcloud blob using number of articles without multiplying by occurences
def normalized_word_blob_by_article_count(rsid_list):
    stats = word_statistics(rsid_list)
    #return " ".join(("{} ".format(noun) * stats[noun][1] for noun in stats))
    blob = ""
    for noun in stats:
        blob += "{} ".format(noun) * stats[noun][1] # noun * article_count
    return blob
    #word_counts = {}
    #pmid_list = (
    #    str(article.pmid) for article in session.query(Article).join(
    #        RSIDCitation,
    #        Article.pmid == RSIDCitation.pmid
    #    ).filter(
    #        RSIDCitation.rsid.in_(rsid_list)
    #    ).distinct().all()
    #)
    #for pmid in pmid_list:
    #    unique_nouns = set(article_noun_mapping[pmid])
    #    for noun in unique_nouns:
    #        if noun in word_counts:
    #            word_counts[noun] += 1
    #        else:
    #            word_counts[noun] = 1
    #normalized_blob = ""
    #for word in word_counts:
    #    normalized_blob += "{} ".format(word) * word_counts[word]
    #return normalized_blob


# Generate a wordcloud png file from a list of rsids
def create_wordcloud_from_rsids(rsid_list, weights, output_path=None):#normalization_type, output_path=None):
    #full_text = {
    #    "default": word_blob,
    #    "rsid": normalized_word_blob_by_rsid,
    #    "article": normalized_word_blob_by_articles,
    #    "article_count": normalized_word_blob_by_article_count,
    #}[normalization_type](rsid_list)

    #full_text = weighted_word_blob(rsid_list, weights)
    stats = word_statistics(rsid_list, weights)
    full_text = word_blob(stats)

    word_cloud = WordCloud(width=1600, height=800).generate(full_text)
    plt.figure(figsize=(20, 10), facecolor='k')
    plt.imshow(word_cloud)
    plt.axis("off")
    plt.tight_layout(pad=0)
    plt.show()
    plt.savefig(output_path, facecolor='k', bbox_inches='tight')

    return stats


# Function results in IndexError if no results are found
def generate_wordcloud(rsid_list, weights):#normalization_type="default"):
    rsid_string = " ".join(rsid_list) + "".join(weights)#normalization_type
    hash_object = hashlib.sha1(rsid_string.encode())
    hashed_file_name = "{}.png".format(hash_object.hexdigest())
    hashed_file_path = os.path.join(WORDCLOUD_STORAGE, hashed_file_name)
    try:
        stats = create_wordcloud_from_rsids(rsid_list, weights, hashed_file_path)#normalization_type, hashed_file_path)
        return hashed_file_name, stats
    except IndexError:
        return False
