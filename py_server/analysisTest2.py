# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 10:41:47 2014
1) app.after_request necesario para implementar política similar a CORS
2) Todos los métodos deben encontrarse entre Flask(__name__) y app.run()

@author: rodri
"""
# --------------------- LIBRARIES -----------------
from flask import Flask, jsonify, request
from flask import redirect, url_for #for uploading files
from werkzeug.utils import secure_filename
import os
import time
import string

#Our methods
import suffixSearch as ss
import motifSearch as ms
import annotations as ann
import helpers
import re

import pickle #i tried cPickle but is way slower!
    

#----------------------- UPLOADS -----------------------
#%%from http://flask.pocoo.org/docs/patterns/fileuploads/
UPLOAD_FOLDER = './genomes' #wherever we run analysis.py
ALLOWED_EXTENSIONS = set(['txt', 'wig', 'bw'])

app=Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
#%%
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS
           
           
@app.route("/")
def hello_world():
	return "Hello, world tangerine!"