#!/usr/bin/env python3
from flask import Flask, request, render_template
from visualization import word_blob
application = Flask(__name__)


@application.route("/")
def index():
    return render_template('index.html')


@application.route("/results/", methods=["GET", "POST"])
def results():
    rsid_string = request.form["rsids"].strip("rs")
    rsid_list = rsid_string.split()
    wordcloud_text = word_blob(rsid_list)
    return render_template('wordcloud.html', wordcloud_text=wordcloud_text)


@application.route("/test")
def test():
    test_rsids = [328]
    blob = word_blob(test_rsids)
    return render_template('wordcloud.html', wordcloud_text=blob)

if __name__ == "main":
    application.run()
