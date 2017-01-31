#!/usr/bin/env python3
from flask import Flask
application = Flask(__name__)


@application.route("/")
def index():
    return "PhenVar!"

if __name__ == "main":
    application.run()
