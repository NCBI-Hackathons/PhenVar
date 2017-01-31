#!/usr/bin/env python3
from flask import Flask, render_template
application = Flask(__name__)


@application.route("/")
def index():
    return render_template('index.html', name=index)

if __name__ == "main":
    application.run()
