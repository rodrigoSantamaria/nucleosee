# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 10:41:47 2014
1) app.after_request necesario para implementar política similar a CORS
2) Todos los métodos deben encontrarse entre Flask(__name__) y app.run()

@author: rodri
"""
# --------------------- INTERNAL METHODS -----------------
#%% READ WIG
def readWig(path="/Users/rodri/Documents/investigacion/IBFG/nucleosomas/Mei3h_center.wig"):
    import numpy
    import time
    t0=time.clock()
    f=open(path)
    seq=f.readlines()
    print '{} s in reading'.format(time.clock()-t0) #5 secs
    t0=time.clock()
    chsize=[]
    cont=0
    for i in range(len(seq)):
        s=seq[i]
        if s[0]=='t' and i>0:#new chromosome
         chsize.append(cont-1)
         cont=0
        else:
            cont=cont+1
    chsize.append(cont-1)
    print '{} s in computing sizes'.format(time.clock()-t0) #39 seqs, go to numpy.array
    t0=time.clock()
    ch=[]
    cont=0
    for i in chsize:
        print i
        cont=cont+2
        chi=numpy.empty(i,dtype=float)
        for j in range(0,i-1):
            chi[j]=round(float(seq[cont+j]),2)
        cont=cont+i
        ch.append(chi)
    print '{} s in formatting'.format(time.clock()-t0) #39 seqs, go to numpy.array
    return ch

#%% DISCRETIZATION
"""Given a numerical sequence seq, this method binarizes based on the average 
and standard deviations on windows of size windowSize. Binarzation is done
in categories a to z (z the larger), as many as detailed by numBins"""
def discretize(seq, windowSize,numBins=5):
    import numpy as np
    alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
    alphabet=alphabetTotal[:numBins]   
    dseq=[]
    sm=np.mean(seq)
    ssd=np.std(seq)
    for i in range(0, len(seq),windowSize):
        im=np.mean(seq[i:i+windowSize])
        #print max(0,min(len(alphabet)/2+int(np.round((im-sm)/ssd)),numBins))
        dseq.append(alphabet[max(0,min(len(alphabet)/2+int(np.round((im-sm)/ssd)),numBins-1))])
    return dseq
    
# --------------------------------------------
#%%
from flask import Flask, jsonify, request
from flask import redirect, url_for #for uploading files
from werkzeug.utils import secure_filename
import os


#%%from http://flask.pocoo.org/docs/patterns/fileuploads/
UPLOAD_FOLDER = '/Users/jonatan/WebStormProjects/seqview/py/genomes' #maybe an absolute path??
#UPLOAD_FOLDER = '.' #wherever we run analysis.py
ALLOWED_EXTENSIONS = set(['txt', 'wig'])

app=Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS
           
           
@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        file = request.files['file']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            print filename
            if filename in os.listdir(os.path.join(app.config['UPLOAD_FOLDER'])):
                print '{} already uploaded'.format(filename)    #TODO: by now we avoid resubmissions (for tests)
            else:
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            #return redirect(url_for('uploaded_file', filename=filename))
            return jsonify(path=os.path.join(app.config['UPLOAD_FOLDER'], filename))
#%%
from flask import send_from_directory

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)
      
@app.route('/testUpload', methods=['GET'])
def testUpload():
    filename=""+request.args.get("filename")
    if filename in os.listdir(os.path.join(app.config['UPLOAD_FOLDER'])):
        return jsonify(response=os.path.join(app.config['UPLOAD_FOLDER'], filename))
    else:
        return jsonify(response='not found')


@app.route('/testUpload2', methods=['GET'])
def testUpload2():
    filename=""+request.args.get("filename")
    if filename in os.listdir(os.path.join(app.config['UPLOAD_FOLDER'])):
        return jsonify(response=os.path.join(app.config['UPLOAD_FOLDER'], filename))
    else:
        return jsonify(response='not found')
      
@app.route("/test")
def test():
    default=""+request.args.get("default") #the only way I found to do it... (no default value can be set)
    default=default.split(",")
    return jsonify(result=len(default))

    
#@app.route("/stats")
#def stats(seq=[0,0,0]):
#    import numpy as np
#    seq=""+request.args.get("seq")
#    seq= [float(x) for x in seq.split(",")]
#    return jsonify(mean=np.mean(seq), sdev=np.std(seq))
#   
#
#@app.route("/discretize")
#def discretize(seq=[0,0,0], windowSize=2,numBins=5):
#    import numpy as np
#    alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
#    alphabet=alphabetTotal[0:numBins]   
#    windowSize=int(request.args.get("windowSize"))
#    numBins=int(request.args.get("numBins"))
#    seq=""+request.args.get("seq")
#    seq= [float(x) for x in seq.split(",")]
#    dseq=[]
#    sm=np.mean(seq)
#    ssd=np.std(seq)
#    for i in range(0, len(seq),windowSize):
#        im=np.mean(seq[i:i+windowSize])
#        dseq.append(alphabet[max(0,min(len(alphabet)/2+int(np.round((im-sm)/ssd)),numBins-1))])
#    print len(dseq)
#    return jsonify(result=dseq)
#
#@app.route("/normalize")
#def normalize(path="/Users/rodri/WebstormProjects/untitled/py/genomes/dwtMini2.wig"):
#    import numpy as np
#    import string
#    path=""+request.args.get("path")
#    f=open(path)
#    seq=f.readlines()
#    del seq[0:2]
#    seq=[float(string.replace(x, "\n", "")) for x in seq]
#    m=np.mean(seq)
#    sd=np.std(seq)
#    nseq=[(x-m)/sd for x in seq]
#    return jsonify(result=nseq, maximum=np.max(nseq), minimum=np.min(nseq), mean=np.mean(seq), sdev=np.std(seq))


"""
Preprocesa un fichero wig para su visualización
path       ruta a los datos (formato wig, de momento sólo se toma el primer track)
windowSize tamaño de la ventana para hacer la discretización (def 100)
numBins    númer de valores discretos para la discretizacion (def 5)
retorna    objeto JSON con los siguientes campos:
               result   valores normalizados (sin discretizar) sampleados hasta maxSize
               fullLength  longitud total de los datos iniciales/normaliados
               maximum,minimum,mean,sdev   medidas estadísticas de los datos normalizados
               dseq        datos discretizados por la función discretize()
"""
@app.route("/preprocess")
def preprocess(path="/Users/rodri/WebstormProjects/untitled/py/genomes/dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000):
    import numpy as np

    #0) read
    print 'reading...'
    seq=readWig(str(request.args.get("path")))
    seq=seq[0] #(TODO: by now, only first chromosome)
    #1) normalize 
    print 'normalizing...'
    m=np.mean(seq)
    sd=np.std(seq)
    nseq=[(x-m)/sd for x in seq]
    #2) discretize
    print 'discretizing...'
    windowSize=int(request.args.get("windowSize"))
    numBins=int(request.args.get("numBins"))
    maxSize=int(request.args.get("maxSize"))
    dseq=discretize(seq, windowSize, numBins)
    print 'done!'
    #TODO: broken pipe returning the object if whole 
            #Although it is not always happening, it's clearly not sensible
            #to send the full sequence. We sample instead 100K values 
    #return jsonify(result=nseq[:min(70000,len(nseq))], fullLength=len(nseq), maximum=np.max(nseq), minimum=np.min(nseq), mean=m, sdev=sd, dseq=dseq)
    #return jsonify(result=[nseq[x] for x in range(0,len(seq),max(1,len(seq)/maxSize))], fullLength=len(nseq), maximum=np.max(nseq), minimum=np.min(nseq), mean=m, sdev=sd, dseq=dseq)
    return jsonify(result=[seq[x] for x in range(0,len(seq),max(1,len(seq)/maxSize))], fullLength=len(seq), maximum=np.max(seq), minimum=np.min(seq), mean=m, sdev=sd, dseq=dseq)

#def preprocess(path="/Users/rodri/WebstormProjects/untitled/py/genomes/dwtMini2.wig", windowSize=500, numBins=5):
#    import numpy as np
#    import string
#    #0) read
#    path=""+request.args.get("path")
#    f=open(path)
#    seq=f.readlines()
#    del seq[0:2]
#    #1) normalize
#    seq=[float(string.replace(x, "\n", "")) for x in seq]
#    m=np.mean(seq)
#    sd=np.std(seq)
#    nseq=[(x-m)/sd for x in seq]
#    #2) discretize
#    windowSize=int(request.args.get("windowSize"))
#    numBins=int(request.args.get("numBins"))
#    alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
#    alphabet=alphabetTotal[0:numBins]   
#    dseq=[]
#    sm=np.mean(nseq)
#    ssd=np.std(nseq)
#    for i in range(0, len(nseq),windowSize):
#        im=np.mean(nseq[i:i+windowSize])
#        dseq.append(alphabet[max(0,min(len(alphabet)/2+int(np.round((im-sm)/ssd)),numBins-1))])
#    return jsonify(result=nseq, maximum=np.max(nseq), minimum=np.min(nseq), mean=sm, sdev=ssd, dseq=dseq)

#%%from http://mortoray.com/2014/04/09/allowing-unlimited-access-with-cors/
@app.after_request
def add_cors(resp):
    """ Ensure all responses have the CORS headers. This ensures any failures 
        are also accessible by the client. """
    resp.headers['Access-Control-Allow-Origin'] = request.headers.get('Origin')
    resp.headers['Access-Control-Allow-Credentials'] = 'true'
    resp.headers['Access-Control-Allow-Methods'] = 'POST, OPTIONS, GET'
    resp.headers['Access-Control-Allow-Headers'] =  request.headers.get('Access-Control-Request-Headers', 'Authorization' )
    # set low for debugging
    if app.debug:
        resp.headers['Access-Control-Max-Age'] = '1'
    return resp
    
if __name__ == '__main__':
    app.run(debug=True) 
    
    
    