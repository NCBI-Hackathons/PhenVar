from lanpros import *
from wordcloud import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def create_wordcloud(normalized_all_counts, outfig):

    word_cloud_list = [(key + ' ') * int(round(normalized_all_counts[key], 4) * 10000)
                       for key in normalized_all_counts.keys()]
    word_cloud_text = '-'.join(word_cloud_list)
    word_cloud = WordCloud().generate(word_cloud_text)
    fig = plt.imshow(word_cloud)
    plt.axis("off")
    plt.savefig(outfig)
