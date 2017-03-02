#!/usr/bin/env python3
from flask import Flask, request, render_template
from datetime import datetime
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
    # Log user's IP address
    with open("logs/visits.log", "a") as logfile:
        log_string = "{}\t{}\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), request.environ["REMOTE_ADDR"])
        logfile.write(log_string)
    rsid_string = request.form["rsids"].strip("rs")
    rsid_list = rsid_string.split()
    normalization_type = request.form["normalization_type"]
    wordcloud_file_name = generate_wordcloud(rsid_list, normalization_type)
    return render_template('results.html', wordcloud_file_name=wordcloud_file_name)

@application.route("/about/")
def about():
    about_markdown = open('markdown/about.md', 'r').read()
    return render_template('about.html', about_markdown=about_markdown)

@application.route("/acknowledgements/")
def acknowledgements():
    acknowledgements_markdown = open('markdown/acknowledgements.md', 'r').read()
    return render_template('acknowledgements.html', acknowledgements_markdown=acknowledgements_markdown)

if __name__ == "main":
    application.run()

