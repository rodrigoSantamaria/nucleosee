# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 10:41:47 2014

@author: rodri
"""
# --------------------- LIBRARIES -----------------
from flask import Flask, jsonify, request
from flask import redirect, url_for #for uploading files
from werkzeug.utils import secure_filename
import os
import time
import string    

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
           
           
@app.route("/test")
def hello_world():
	return "Hello, world tangerine!"
