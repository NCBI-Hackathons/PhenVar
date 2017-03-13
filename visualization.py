from wordcloud import *
from database import session, RSIDCitation, Article, article_noun_mapping
from settings import WORDCLOUD_STORAGE
import matplotlib.pyplot as plt
import hashlib
import os


# Create table - NOUN|RSID's|#Articles
def word_statistics(rsid_list, weights):
    # dictionary in form {"noun": [word_count, article_count, rsids, strength]}
    word_stats = {}
    results = session.query(Article, RSIDCitation).join(
        RSIDCitation,
        Article.pmid == RSIDCitation.pmid
    ).filter(
        RSIDCitation.rsid.in_(rsid_list)
    ).all()
    recorded_pmids = []
    if "rsid_balance" in weights:
        for rsid in rsid_list:
            rsid_stats = {}
            pmid_list = list((str(result.Article.pmid) for result in results if result.RSIDCitation.rsid == int(rsid)))
            normalization_factor = 1000/len(pmid_list)
            for pmid in pmid_list:
                nouns = article_noun_mapping[pmid]
                for noun in set(nouns):
                    if noun not in rsid_stats:
                        rsid_stats[noun] = [0, 0, [], 0]
                    rsid_stats[noun][0] += nouns.count(noun)
                    rsid_stats[noun][1] += 1
            for noun in rsid_stats:
                noun_stats = rsid_stats[noun]
                noun_strength = normalization_factor
                weights_dictionary = {
                    "word_count": noun_stats[0],
                    "article_count": noun_stats[1],
                    "rsid_count": 1,
                    "rsid_balance": 1,
                }
                for weight_label in weights:
                    noun_strength *= weights_dictionary[weight_label]
                noun_stats[3] = int(noun_strength)
                if noun not in word_stats:
                    word_stats[noun] = noun_stats
                else:
                    word_stats[noun][3] += int(noun_strength)
                word_stats[noun][2].append(rsid)
        # Additional weighting by rsid count will require processing after all rsid's have been processed
        if "rsid_count" in weights:
            for noun in word_stats:
                word_stats[noun][3] *= len(word_stats[noun][2])
    else:
        for result in results:
            pmid = str(result.Article.pmid)
            rsid = result.RSIDCitation.rsid
            nouns = article_noun_mapping[pmid]
            counting = pmid not in recorded_pmids
            for noun in set(nouns):
                if counting:
                    if noun not in word_stats:
                        word_stats[noun] = [0, 0, [], 1]
                    word_stats[noun][0] += nouns.count(noun)
                    word_stats[noun][1] += 1
                if rsid not in word_stats[noun][2]:
                    word_stats[noun][2].append(rsid)
            if counting:
                recorded_pmids.append(pmid)

        for noun in word_stats:
            noun_stats = word_stats[noun]
            weights_dictionary = {
                "word_count": noun_stats[0],
                "article_count": noun_stats[1],
                "rsid_count": len(noun_stats[2]),
            }
            for weight_label in weights:
                word_stats[noun][3] *= weights_dictionary[weight_label]
    return word_stats


def associated_articles_rsids(noun):
    results = {}
    pmid_list = [pmid for pmid in article_noun_mapping if noun in article_noun_mapping[pmid]]
    for pmid in pmid_list:
        results[pmid] = list([
            str(result[0]) for result in session.query(RSIDCitation.rsid).filter_by(pmid=int(pmid)).distinct().all()
        ])
    return results


# Generate a wordcloud png file from a list of rsids
def create_wordcloud_from_stats(stats, output_path=None):#rsid_list, weights, output_path=None):#normalization_type, output_path=None):
    #stats = word_statistics(rsid_list, weights)
    frequencies = []
    for noun in stats:
        frequencies.append((noun, stats[noun][3]))
    word_cloud = WordCloud(
        width=1600,
        height=800,
        #relative_scaling=0.75,
    ).generate_from_frequencies(frequencies=frequencies)
    plt.figure(figsize=(20, 10), facecolor='k')
    plt.imshow(word_cloud)
    plt.axis("off")
    plt.tight_layout(pad=0)
    plt.show()
    plt.savefig(output_path, facecolor='k', bbox_inches='tight')
    #return stats


# Function results in IndexError if no results are found
def generate_wordcloud(stats, rsid_list, weights): #rsid_list, weights):
    rsid_string = " ".join(rsid_list) + "".join(weights)
    hash_object = hashlib.sha1(rsid_string.encode())
    hashed_file_name = "{}.png".format(hash_object.hexdigest())
    hashed_file_path = os.path.join(WORDCLOUD_STORAGE, hashed_file_name)
    try:
        create_wordcloud_from_stats(stats, hashed_file_path)#rsid_list, weights, hashed_file_path)
        return hashed_file_name
    except IndexError:
        return False
