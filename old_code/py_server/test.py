# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 10:53:29 2014

@author: rodri
"""

from flask import Flask, jsonify, request, make_response

app=Flask(__name__)



@app.route("/index")
def index():
    return jsonify(result=7)


if __name__ == '__main__':
    app.run(debug=True) 
