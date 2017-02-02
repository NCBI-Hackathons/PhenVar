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
    wordcloud_file_name = generate_wordcloud(rsid_list)
    return render_template('results.html', wordcloud_file_name=wordcloud_file_name)


if __name__ == "main":
    application.run()
