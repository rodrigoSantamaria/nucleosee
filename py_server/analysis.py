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
#sys.path.append("/Users/rodri/WebstormProjects/seqview/py_server")
import suffixSearch as ss
import annotations as ann
import helpers as h


# --------------------- INTERNAL METHODS -----------------
#%% -----------  READ WIG --------------
def readWig(path="/Users/rodri/Documents/investigacion/IBFG/nucleosomas/Mei3h_center.wig"):
    import numpy
    import re
    import time
    t0=time.clock()
    f=open(path)
    seq=f.readlines()
    print '{} s in reading'.format(time.clock()-t0) #2 secs
    t0=time.clock()
    chsize=[]
    cont=0
    for i in xrange(len(seq)):
        s=seq[i]
        if s[0]=='t' and i>0:#new chromosome
         chsize.append(cont-1)
         cont=0
        else:
            cont=cont+1
    chsize.append(cont-1)
    print '{} s in computing sizes'.format(time.clock()-t0) #6 seqs, go to numpy.array
    t0=time.clock()
    ch={}
    cont=0
    name=re.sub("\n", "", re.sub(" .*$", "", re.sub("^.*chrom=", "", seq[cont+1])))
    ch[name]=[]
    print "name is ", name
    for i in chsize:
        print i
        cont=cont+2
        chi=numpy.empty(i,dtype=float)
        for j in xrange(0,i-1):
            #chi[j]=round(float(seq[cont+j]),2)
            chi[j]=seq[cont+j]
        ch[name].append(chi)
        cont=cont+i
        if(cont<len(seq)):
            name=re.sub("\n", "",re.sub(" .*$", "", re.sub("^.*chrom=", "", seq[cont])))
            print "name is ", name
            ch[name]=[]
    print '{} s in formatting'.format(time.clock()-t0) #30 seqs, go to numpy.array
    return ch
    
#tal=readWig("/Users/rodri/WebstormProjects/seqview/py_server/genomes/jpiriz/23479_h90_wlt_mean.wig")

#%% ------------------ DISCRETIZATION -------------------
"""Given a numerical sequence seq, this method binarizes based on the average 
and standard deviations on windows of size windowSize. Binarzation is done
in categories a to z (z the larger), as many as detailed by numBins"""
def discretize(seq, windowSize, minimo, maximo, numBins=5):
    import numpy as np
    alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
    alphabet=alphabetTotal[:numBins]   
    dseq=[]
    factor=(numBins-1)/(maximo-minimo)
    sseq=np.split(np.array(seq[:windowSize*(len(seq)/windowSize)]), len(seq)/windowSize)
    mseq=np.mean(sseq, axis=1, keepdims=True)
    for im in mseq:
        dseq.append(alphabet[(int)(factor*(im-minimo))])
    return dseq 

#def discretize(seq, windowSize, minimo, maximo, numBins=5):
#    import numpy as np
#    alphabetTotal=['a','b','c','d','e', 'f', 'g','h','i','j','k','l','m','n','o','p','q','r','s','t']
#    alphabet=alphabetTotal[:numBins]   
#    factor=(numBins-1)/(maximo-minimo)
#    sseq=np.split(np.array(seq[:windowSize*(len(seq)/windowSize)]), len(seq)/windowSize)
#    mseq=np.mean(sseq, axis=1, keepdims=True)
#    def disc(a):
#        return alphabet[(int)(factor*(a-minimo))]
#    return list(np.apply_along_axis(disc, 1,mseq))
 
#import time
#t0=time.clock()   
#tal=discretize(seq,30, minimum, maximum)
#print '{}'.format((time.clock()-t0))

#----------------------- UPLOADS -----------------------
#%%from http://flask.pocoo.org/docs/patterns/fileuploads/
#UPLOAD_FOLDER = '/Users/rodri/WebstormProjects/seqview/py_server/genomes' #maybe an absolute path??
UPLOAD_FOLDER = './genomes' #wherever we run analysis.py
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
def preprocess(filename="dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000, track=None):
    import numpy as np
    import time
    global data
    global session

    t00=time.clock()
    #0) read
    print 'reading...'
    path=os.path.join(app.config['UPLOAD_FOLDER'],user,str(request.args.get("filename")))
    track=request.args.get("track")
    if track=="None":
        genome=readWig(path)
        track=genome.keys()[0]
        seq=genome[track][0]
    else:
        genome=data["genome"]
        seq=seq[track]

    #1) normalize 
    t0=time.clock()
    print 'computing statistics...'
    m=np.mean(seq)
    sd=np.std(seq)
    
    upperlim=m+3*sd#avoid outliers? testing
    seq=np.clip(seq,0,upperlim)

    maximum=np.max(seq)
    minimum=np.min(seq)
    print 'max and min in {}s'.format((time.clock()-t0))


    #2) discretize
    t0=time.clock()
    print 'discretizing...'
    windowSize=int(request.args.get("windowSize"))
    numBins=int(request.args.get("numBins"))
    maxSize=int(request.args.get("maxSize"))
    dseq=discretize(seq, windowSize, minimum, maximum, numBins)
    print 'done! in {}s'.format((time.clock()-t0))

    t0=time.clock()
    print 'bwt...'
    t=ss.bwt(''.join(dseq)+"$")
    print 'done! in {}s'.format((time.clock()-t0))

    t0=time.clock()
    print 'sampling...'
    res=[round(seq[x],2) for x in xrange(0,len(seq),max(1,len(seq)/maxSize))]
    #step=max(1,len(seq)/maxSize)
    #res=[np.mean(seq[x:x+step]) for x in range(0,len(seq),step)] #TODO: round?
    print 'done! in {}s'.format((time.clock()-t0))
    print 'sample length {} and an element {}'.format(len(res), res[44])

    t0=time.clock()
    print 'loading annotations...'
    print(h.gffPath(ch=track))
    dataGFF=ann.gff(h.gffPath(ch=track))
    print 'time in GFF:',(time.clock()-t0),'s'
    t0=time.clock()
    dataGO=ann.go()
    print 'time in GO:',(time.clock()-t0),'s'
    t0=time.clock()
    dataGOA=ann.goa()
    print 'time in GOA:',(time.clock()-t0),'s'
    #dataFASTA=fasta(1)
    print "done! ... annotations takes {}".format((time.clock()-t0))
    data={"genome": genome, "seq":seq, "res":res, "fullLength":len(seq), "maximum":maximum, "minimum":minimum,
          "mean":m, "stdev":sd, "dseq":dseq, "bwt":t, "gff":dataGFF,
          "go":dataGO, "goa":dataGOA}
    session[user]=data
    print 'whole preprocess takes {}s'.format((time.clock()-t00))
    return jsonify(seq=res, fullLength=len(seq), maximum=maximum, minimum=minimum, mean=m, stdev=sd, dseq=dseq, chromosomes=genome.keys())


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
    print "{} before".format(pattern)
    pattern=h.convertString(pattern)
    print "{} after".format(pattern)
#    t0=time.clock()
#    data=loadSession(session)
#    print "Session load takes {}".format((time.clock()-t0))
    t=data["bwt"]
    print t["firstOccurrence"]
    if(False in [x in t["firstOccurrence"] for x in set(pattern)]):
        return jsonify(response="There are characters in pattern that do not correspond to the sequence characters: {}".format(t["firstOccurrence"].keys()))
    else:
        t0=time.clock()
        match=ss.bwMatchingV8("".join(data["dseq"]), pattern, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"],1000, d)
        print "Search takes {}".format((time.clock()-t0))
        return jsonify(points=(str)(match), sizePattern=len(pattern))



#%%
@app.route("/getPartSeq")
def getPartSeq(start=0, end=0):
    global data

    start=int(request.args.get("start"))
    end=int(request.args.get("end"))
    seq=data["seq"]
    part=seq[start:end]
    return jsonify(partSeq=part)



#%%
@app.route("/annotations")
def annotations(positions=[], window=1000, types=["any"]):
    import time
    t0=time.clock()
    global data
    window=int(request.args.get("window"))
    print window
    pos=eval(request.args.get("positions"))
    print pos
    types=eval(request.args.get("types"))
    print types
    res=ann.annotate(pos, data["gff"], types, window)
    print "Annotations take {}".format((time.clock()-t0))
    return jsonify(response=res)

#%%
@app.route("/annotationsGOA")
#annotations is the result of calling annotations
#types by now is any. Maybe in the future we divide in BP,MF,CC
def annotationsGOA(annotations={}, types=["any"]):
    global data
    
    a=eval(request.args.get("annotations"))
    print 'annotation size: {}'.format(len(a))
    agon=ann.annotateGOnames(a, data["goa"], data["go"])
    return jsonify(response=agon)
    
#%%
@app.route("/enrichmentGO")
#annotations is the result of calling annotations
#correction is for multiple hipothesis and can be 'none', 'bonferroni', 'fdr' or 'fwer'
#alpha is the threshold, applied as is with 'none' correction or with the specified one
def enrichmentGO(annotations={}, correction="none", alpha=0.01):
    global data
    
    annotations=eval(request.args.get("annotations"))
    alpha=float(request.args.get("alpha"))
    correction=request.args.get("correction")
    
    gis=set()
    for x in annotations.keys():
        for y in annotations[x]:
            gis.add(y["id"])
    ego={}
    ego=ann.enrichmentFisher(gis, data["goa"], alpha, correction)
    for k in ego.keys():
        if(data["go"].has_key(k) and data["go"][k]!="undefined"):
            ego[k]["go_name"]=data["go"][k]
    return jsonify(response=ego)


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
    
    
    