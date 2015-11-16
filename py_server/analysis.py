# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 10:41:47 2014
1) app.after_request necesario para implementar política similar a CORS
2) Todos los métodos deben encontrarse entre Flask(__name__) y app.run()

@author: rodri
"""
"""PRUEBA GIT"""
# --------------------- LIBRARIES -----------------
from flask import Flask, jsonify, request
from flask import redirect, url_for #for uploading files
from werkzeug.utils import secure_filename
import os
import time

#BWT searches
import sys 
sys.path.append("/Users/rodri/WebstormProjects/seqview/py_server")
import suffixSearch as ss


# --------------------- INTERNAL METHODS -----------------
#%% -----------  READ WIG --------------
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

#%% ------------------ DISCRETIZATION -------------------
"""Given a numerical sequence seq, this method binarizes based on the average 
and standard deviations on windows of size windowSize. Binarzation is done
in categories a to z (z the larger), as many as detailed by numBins"""
def discretize(seq, windowSize,numBins=5):
    import numpy as np
    alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
    alphabet=alphabetTotal[:numBins]   
    dseq=[]
    #sm=np.mean(seq)
    #ssd=np.std(seq)
    maximo=max(seq)
    minimo=min(seq)
    for i in range(0, len(seq),windowSize):
        im=np.mean(seq[i:i+windowSize])
        dseq.append(alphabet[(int)(np.round((numBins-1)*(im-minimo)/(maximo-minimo)))])
        #dseq.append(alphabet[max(0,min(len(alphabet)/2+int(np.round((im-sm)/ssd)),numBins-1))])
    return dseq
    
# --------------------------------------------


#----------------------- UPLOADS -----------------------
#%%from http://flask.pocoo.org/docs/patterns/fileuploads/
UPLOAD_FOLDER = '/Users/rodri/WebstormProjects/seqview/py_server/genomes' #maybe an absolute path??
#UPLOAD_FOLDER = '.' #wherever we run analysis.py
ALLOWED_EXTENSIONS = set(['txt', 'wig'])

app=Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS
           
           
@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    global user
    if request.method == 'POST':
        file = request.files['file']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            print filename
            root=os.path.join(app.config['UPLOAD_FOLDER'])
            if (user in os.listdir(root))==False:                 #create directory for this user
                os.mkdir(os.path.join(root,user))
            if filename in os.listdir(os.path.join(root,user)):
                print '{} already uploaded'.format(filename)    #TODO: by now we avoid resubmissions (for tests)
            else:
                print 'uploading...'
                file.save(os.path.join(root, user, filename))
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
    cpath=os.path.join(app.config['UPLOAD_FOLDER'],user)
    if os.path.exists(cpath) and filename in os.listdir(cpath):
        import hashlib
        codeEx=""+request.args.get("md5")
        codeIn=hashlib.md5(open(os.path.join(cpath,filename), 'rb').read()).hexdigest()
        if(codeEx==codeIn):
            return jsonify(response="file exists")
        else:
            return jsonify(response='outdated version')
    else:
        return jsonify(response='not found')
    
      
#@app.route("/test")
#def test():
#    default=""+request.args.get("default") #the only way I found to do it... (no default value can be set)
#    default=default.split(",")
#    return jsonify(result=len(default))

#%%------------- SESSION MANAGEMENT
#data must be a dic with dseq and bwt
def saveSession(path="/Users/rodri/WebstormProjects/untitled/py_server/genomes/dwtMini2.pkl", data=""):
    import cPickle as pickle
    f=open(path, "wb")
    pickle.dump(data, f)
    return 0
    
#data must be a dic with dseq and bwt TODO: check times here
def loadSession(path="/Users/rodri/WebstormProjects/untitled/py_server/genomes/dwtMini2.pkl"):
    import cPickle as pickle
    f=open(path, "r")
    session=pickle.load(f)
    return session
#%%
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
def preprocess(filename="dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000, track=0):
    import numpy as np
    global data
    global session
    #0) read
    print 'reading...'
    path=os.path.join(app.config['UPLOAD_FOLDER'],user,str(request.args.get("filename")))
    track=int(request.args.get("track"))
    seq=readWig(path)
    seq=seq[track] #(TODO: by now, only first chromosome)
    #1) normalize 
    print 'computing statistics...'
    m=np.mean(seq)
    sd=np.std(seq)
    maximum=np.max(seq)
    minimum=np.min(seq)
    #nseq=[(x-m)/sd for x in seq]
    #2) discretize
    print 'discretizing...'
    windowSize=int(request.args.get("windowSize"))
    numBins=int(request.args.get("numBins"))
    maxSize=int(request.args.get("maxSize"))
    dseq=discretize(seq, windowSize, numBins)
    print 'done!'
    print 'bwt...'
    t=ss.bwt(''.join(dseq)+"$")
    print 'done!'
    res=[round(seq[x],2) for x in range(0,len(seq),max(1,len(seq)/maxSize))]
    data={"seq":res, "fullLength":len(seq), "maximum":maximum, "minimum":minimum, "mean":m, "stdev":sd, "dseq":dseq, "bwt":t}
    session[user]=data
    return jsonify(seq=res, fullLength=len(seq), maximum=maximum, minimum=minimum, mean=m, stdev=sd, dseq=dseq)


#%% -------------- SEARCHES -----------

"""
Busca un determinado texto con una búsqueda en el array de sufijos creado
De momento estamos usando variables globales (!) para el texto y la estructura BWT
pattern    patrón de búsqueda
d          nº de mutaciones permitidas
retorna    las posiciones dentro de dseq donde aparece el patrón
"""
#TODO: This works, but the session load is costly. Should be improved or substituted by other method
@app.route("/search")
def search(pattern="", d=0):
    global data
    
    d=int(request.args.get("d"))
    pattern=str(request.args.get("pattern"))
#    t0=time.clock()
#    data=loadSession(session)
#    print "Session load takes {}".format((time.clock()-t0))
    t=data["bwt"]
    print t["firstOccurrence"]
    if(False in [x in t["firstOccurrence"] for x in set(pattern)]):
        return "There are characters in pattern that do not correspond to the sequence characters: {}".format(t["firstOccurrence"].keys())
    else:
        t0=time.clock()
        match=ss.bwMatchingV8("".join(data["dseq"]), pattern, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"],1000, d)
        print "Search takes {}".format((time.clock()-t0))
        return (str)(match)#return "hola {} {}".format(pattern, d)



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
    
#%% session info
@app.before_request
def load_passport():
    global data
    global session
    global user
    user=str(request.args.get("user"))
    #password=str(request.args.get("password"))
    #TODO: check password and so on.
    if user in session.keys():
        data=session[user]
        print len(data)
        print data.keys()

#@app.after_request
#def serialize_passport(response):
#    if hasattr(g, "passport"):
#        session["passport_id"] = g.passport.id
#    return response
    

#-------------------- LAUNCH -----------------
if __name__ == '__main__':
    global session
    global data
    global user
    user=""
    session={}
    session["jpiriz"]={}    
    session["rodri"]={}    
    data={}
    app.run(debug=True) 
    
    
    