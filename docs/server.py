#!/usr/bin/env python

import flask

app = flask.Flask(
    __name__, 
    static_folder='_build/html',
)

@app.route('/')
def home():
    return flask.send_from_directory('_build/html',  'index.html')

@app.route('/<path:path>')
def page(path):
    return flask.send_from_directory('_build/html',  path)

if __name__ == "__main__":
    app.run(debug=True)
