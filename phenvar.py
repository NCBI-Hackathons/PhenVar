#!/usr/bin/env python3
from flask import Flask, request, render_template
from flaskext.markdown import Markdown
from visualization import generate_wordcloud
application = Flask(__name__)
Markdown(application)


@application.route("/")
def index():
    home_markdown = open('markdown/home.md', 'r').read()
    return render_template('index.html', home_markdown=home_markdown)


@application.route("/results/", methods=["GET", "POST"])
def results():
    rsid_string = request.form["rsids"].strip("rs")
    rsid_list = rsid_string.split()
    normalization_type = request.form["normalization_type"]
    wordcloud_file_name = generate_wordcloud(rsid_list, normalization_type)
    return render_template('results.html', wordcloud_file_name=wordcloud_file_name)

@application.route("/about/")
def about():
    about_markdown = open('markdown/about.md', 'r').read()
    return render_template('about.html', about_markdown=about_markdown)

if __name__ == "main":
    application.run()
