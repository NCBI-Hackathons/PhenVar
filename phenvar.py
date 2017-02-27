#!/usr/bin/env python3
from flask import Flask, request, render_template
from datetime import datetime
from visualization import generate_wordcloud
application = Flask(__name__)


@application.route("/")
def index():
    return render_template('index.html')


@application.route("/results/", methods=["GET", "POST"])
def results():
    # Log user's IP address
    with open("logs/visits.log", "a") as logfile:
        log_string = "{}\t{}\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), request.environ["REMOTE_ADDR"])
        logfile.write(log_string)
    rsid_string = request.form["rsids"].strip("rs")
    rsid_list = rsid_string.split()
    wordcloud_file_name = generate_wordcloud(rsid_list)
    return render_template('results.html', wordcloud_file_name=wordcloud_file_name)

@application.route("/about/")
def about():
    return render_template('about.html')

if __name__ == "main":
    application.run()
