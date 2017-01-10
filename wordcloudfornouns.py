from lanpros import *
from wordcloud import *
import matplotlib.pyplot as plt

def create_wordcloud(normalized_all_counts):

    word_cloud_list = [(key + ' ') * int(round(normalized_all_counts[key],4)*10000) for key in normalized_all_counts.keys()]
    word_cloud_text = '-'.join(word_cloud_list)
    word_cloud = WordCloud().generate(word_cloud_text)
    plt.imshow(word_cloud)
    plt.axis("off")
    plt.show()
    return


#test implementation
"""
toquery=get_pmids("rs12255372")

abstracts = get_abstracts(toquery)

t = tokenize_abstracts(abstracts)
q = tagged_abstracts(t)
n = extract_nouns(q)
create_wordcloud(n)
"""
