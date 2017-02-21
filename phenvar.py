#!/usr/bin/env python3
from flask import Flask, request, render_template
from visualization import generate_wordcloud
application = Flask(__name__)


@application.route("/")
def index():
    return render_template('index.html')


@application.route("/results/", methods=["GET", "POST"])
def results():
    rsid_string = request.form["rsids"].strip("rs")
    rsid_list = rsid_string.split()
    normalization_type = request.form["normalization_type"]
    wordcloud_file_name = generate_wordcloud(rsid_list, normalization_type)
    return render_template('results.html', wordcloud_file_name=wordcloud_file_name)

@application.route("/about/")
def about():
    return render_template('about.html')

if __name__ == "main":
    application.run()
