#!/usr/bin/env python

import subprocess
import flask

app = flask.Flask(
    __name__, 
    static_folder='_build/html',
)

@app.route('/')
def home():
    return flask.send_from_directory('_build/html',  'index.html')

@app.route('/sync', methods=['GET', 'POST'])
def sync():
    p = subprocess.run(
        'GIT_DIR=../../git GIT_WORK_TREE=.. git pull origin master',
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    r = p.stdout.decode()
    if r != "Already up to date.\n":
        r += restart()
    return f'<p>{r}</p>'


@app.route('/restart')
def restart():
    root = flask.request.base_url.split('/')[2]
    p = subprocess.run(
        f'sudo supervisorctl restart {root}',
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out = p.stdout.decode()
    err = p.stderr.decode()
    return f"<ul><li>{out}</li> <li>{err}</li></ul>"

@app.route('/<path:path>')
def page(path):
    return flask.send_from_directory('_build/html',  path)

if __name__ == "__main__":
    app.run(debug=True)
