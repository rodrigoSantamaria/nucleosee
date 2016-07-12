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

#Our methods
import suffixSearch as ss
import motifSearch as ms
import annotations as ann
import helpers
import pickle
import re
    

    


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
            print(filename)
            root=os.path.join(app.config['UPLOAD_FOLDER'])
            if (user in os.listdir(root))==False:                 #create directory for this user
                os.mkdir(os.path.join(root,user))
            if filename in os.listdir(os.path.join(root,user)):
                print(filename,'already uploaded')    #TODO: by now we avoid resubmissions (for tests)
            else:
                print('uploading...')
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
    filename=re.sub(r"\..*$", ".pic", filename)
    print(filename)
    cpath=os.path.join(app.config['UPLOAD_FOLDER'],user)
    if os.path.exists(cpath) and filename in os.listdir(cpath):
        codeEx=request.args.get("forceReload")
        if(codeEx=="True"):
            return jsonify(response='outdated version')
        else:
            return jsonify(response="file exists")
    else:
        return jsonify(response='not found')
 
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
maxSize    maximum length of the returning array (to avoid bandwidth overload). Default 100K
stdev      outlier clipping above this number of standard deviations from the mean. Default 3
track      in the case of several tracks, which one to return. Default 'None' returns the first track on alphabetic order
reload     if the file must be reloaded (True) or not, if a previous preprocessed version exists in pickle format (False). Default False
retorna    objeto JSON con los siguientes campos:
               result   valores normalizados (sin discretizar) sampleados hasta maxSize
               fullLength  longitud total de los datos iniciales/normaliados
               maximum,minimum,mean,sdev   medidas estadísticas de los datos normalizados
               dseq        datos discretizados por la función discretize()
"""
#%%
@app.route("/preprocess")
def preprocess(filename="dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000, stdev=3, track="None", recharge="False"):
    import numpy as np
    global data
    global session

    
    t00=time.clock()
    #0) read
    print('reading...')
    t0=time.clock()
    basePath=os.path.join(app.config['UPLOAD_FOLDER'],user)
    filename=str(request.args.get("filename"))
    path=os.path.join(basePath,filename)
    track=request.args.get("track")
    forceReload=request.args.get("recharge")
    windowSize=int(request.args.get("windowSize"))
    numBins=int(request.args.get("numBins"))
    maxSize=int(request.args.get("maxSize"))
    stdev=float(request.args.get("stdev"))
    picklePath=os.path.join(basePath,re.sub(r"\..*$", ".pic", filename))
    savePickle=False
    
    m={}
    sd={}
    maximum={}
    minimum={}
    dseq={}#discretized (alphanumeric) sequence
    bins={}#thresholds for the binaries
    res={}#reduced sequence
    t={} #burrows-wheeler transform
    dataGFF={}
    dataFASTA={}
    seqd={}
    
    if os.path.isfile(picklePath) and forceReload=="False":
        print("Pickle exists!!!, recharge=", forceReload)
        f=open(picklePath)
        datap=pickle.load(f)
        genome=datap["seq"]
        t=datap["bwt"]
        dseq=datap["dseq"]
        maximum=datap["maximum"]
        minimum=datap["minimum"]
        m=datap["mean"]
        sd=datap["stdev"]
        bins=datap["bins"]
    else:
        genome=helpers.readWig(path)
        savePickle=True


    if track=="None":
        track=sorted(genome.keys())[0]
    
    print("load wig takes",(time.clock()-t0))
   
    #for each separate chromosome: TODO: include ds and bwt into pickle!
    for k in genome.keys():
        tk=time.clock()
        seq=genome[k]
        print("ch", k, "with length", len(seq))
        if(savePickle):
            #1) normalize 
            t0=time.clock()

  
            m[k]=np.mean(seq, dtype=float)
            sd[k]=np.std(seq, dtype=float)
            upperlim=m[k]+stdev*sd[k]#avoid outliers? testing
            seq=np.clip(seq,0,upperlim)
        
            m[k]=np.mean(seq, dtype=float)
            sd[k]=np.std(seq, dtype=float)
            maximum[k]=np.max(seq)
            minimum[k]=np.min(seq)
            print('\\tstats in ',(time.clock()-t0), "s")
        
        
            #2) discretize
            t0=time.clock()
            tmp=helpers.discretize(seq, windowSize, minimum[k], maximum[k], numBins, percentile=True)
            dseq[k]=tmp["dseq"]
            bins[k]=tmp["bins"]
            print('\\tdiscretize in',(time.clock()-t0),' s')
        
            t0=time.clock()
            t[k]=ss.bwt(''.join(dseq[k])+"$")
            print('\tbwt in ',(time.clock()-t0),'s')
    
        t0=time.clock()
        res[k]=list(np.mean(helpers.rolling_window(seq, max(1,len(seq)/maxSize)),-1, dtype=float)) #maybe round?  
        print('\tsampling in',(time.clock()-t0),'s')
    
        t0=time.clock()
        dataGFF[k]=ann.gff(helpers.gffPath(ch=k))
        dataFASTA[k]=ann.fasta(k)
        print('\ttime in annotations (GFF and FASTA):',(time.clock()-t0),'s')
        seqd[k]=seq
        print("processing",k,"takes",(time.clock()-tk))
        
    #3) annotations
    t0=time.clock()
    dataGO=ann.go()
    dataGOA=ann.goa()
    print("done! ... GO annotations takes",(time.clock()-t0))
    
    data={"seq":seqd, "fullLength":len(seq), "maximum":maximum, "minimum":minimum,
          "mean":m, "stdev":sd, "dseq":dseq, "bwt":t, "gff":dataGFF, "res":res,
          "go":dataGO, "goa":dataGOA, "fasta":dataFASTA, "bins":bins, "windowSize":windowSize}
    session[user]=data
    
    if(savePickle):
        tpickle=time.clock()
        print("pickle path:", picklePath)
        datap={"seq":seqd, "dseq":dseq, "bwt":t, "maximum":maximum, "minimum":minimum,
          "mean":m, "stdev":sd, "bins":bins, "windowSize":windowSize}
    
        f=open(picklePath, 'w')
        pickle.dump(datap,f)
        print("serialize data takes",(time.clock()-tpickle))
        #NOTE: should remove the .wig here to avoid double memory space?

    print('whole preprocess takes ',(time.clock()-t00),"s")
    print('returning track',track)
    print('max values are',maximum)
    
    return jsonify(seq=res[track], fullLength=len(genome[track]), maximum=(float)(maximum[track]), minimum=(float)(minimum[track]), mean=(float)(m[track]), stdev=(float)(sd[track]), dseq=dseq[track], bins=bins[track], chromosomes=sorted(genome.keys()))

#%%preprocess(filename="/Users/rodri/Documents/investigacion/IBFG/nucleosomas/dwtMini2.wig")

#%%
@app.route("/getTrack")
def getTrack(track="None"):
    global data

    track=request.args.get("track")
    if track=="None":
        track=data["res"].keys()[0]

    print("returning chromosome",track)
    if("search" in data.keys() == False):
        d["search"]={"points":{}, "sizePattern":0}
    return jsonify(seq=data["res"][track], fullLength=len(data["seq"][track]), maximum=(float)(data["maximum"][track]), minimum=(float)(data["minimum"][track]), mean=(float)(data["mean"][track]), stdev=(float)(data["stdev"][track]), dseq=data["dseq"][track], bins=data["bins"][track], chromosomes=sorted(data["res"].keys()), search=data["search"])

#%% -------------- SEARCHES -----------

"""
Busca un determinado texto con una búsqueda en el array de sufijos creado
De momento estamos usando variables globales (!) para el texto y la estructura BWT
pattern    patrón de búsqueda
d          nº de mutaciones permitidas
geo        si distinto de "none", filtra las búsquedas que no esten en la zona genomica indicada ("gene", "exon", "utr", "ncRNA gene")
retorna    las posiciones dentro de dseq donde aparece el patrón
"""
@app.route("/search")
def search(pattern="", d=0, geo="none", intersect="soft", softMutations="false"):
    global data
    
    print("searching..:")
    t00=time.clock()
    
    d=int(request.args.get("d"))
    pattern=str(request.args.get("pattern"))
    geo=str(request.args.get("geo"))
    intersect=str(request.args.get("intersect"))
    softMutations=str(request.args.get("softMutations"))

    ws=data["windowSize"]
                
    #CASE 1) Numerical range pattern
    interval=helpers.convertRange(pattern)
    if(interval!=-1):
        print("NUMERICAL RANGE")
        points={}
        for k in data["seq"].keys():
            #points[k]=(str)([(int)(interval["start"]/ws)] if (interval["start"]+interval["length"])<len(data["seq"][k]) else [])
            points[k]=(str)([(int)(interval["start"])] if (interval["start"]+interval["length"])<len(data["seq"][k]) else [])
        data["search"]={'points':points, 'sizePattern':((int)(interval["length"]/ws))}
        return jsonify(points=points, sizePattern=((int)(interval["length"]/ws)));
    
    #CASE 2) gene/go name pattern
    t=data["bwt"][data["seq"].keys()[0]]
    patternLetters=t["firstOccurrence"].keys()
    patternLetters.append('+')
    patternLetters.append('*')
    for x in range(1,10):
        patternLetters.append((str)(x)) 
    
    print(patternLetters, "\t",set(pattern))
    if(False in [x in patternLetters for x in set(pattern)]):
        print("GENE OR TERM")       #TODO: GO Term search not implemented yet
        points={}
        sizes={}
        for k in data["seq"].keys():
            oc=ann.searchGene(pattern, data["gff"][k])
            #points[k]=(str)([(int)(y["start"]/ws) for y in oc])
            points[k]=(str)([(int)(y["start"]) for y in oc])
            sizes[k]=(str)([(int)((y["end"]-y["start"])/ws) for y in oc])
        data["search"]={'points':points, 'sizePattern':sizes}
        return jsonify(points=points, sizePattern=sizes);
        
    #CASE 3) bin pattern    
    pattern=helpers.convertString(pattern)
    print("BWT ON ",pattern)
    search={}
    for k in data["seq"].keys():
        if(False in [x in t["firstOccurrence"] for x in set(pattern)]):
            return jsonify(response="error", msg="{}: There are characters in pattern that do not correspond to the sequence characters: {}".format(k, t["firstOccurrence"].keys()))
        else:
            t0=time.clock()
            t=data["bwt"][k]
            search[k]=(ss.bwMatchingV8("".join(data["dseq"][k]), pattern, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"],1000, d))
            print(len("".join(data["dseq"][k])))
            print("Search ",k,"takes",(time.clock()-t0), "and finds",len(search[k]), "occurences")
            if(len(search[k])>10000):
                return jsonify(response="error", msg="Too many occurrences, please narrow your search", points={}, sizePattern=len(pattern))
    print("Search finished in ",(time.clock()-t00))
    
    
    
    
    if(softMutations=="true"):
        ti=time.clock()
        for k in data["seq"].keys():
            search[k]=helpers.filterHard(search[k], data["dseq"][k], pattern)
        print("Soft mutation filtering takes ", (time.clock()-ti))
    
    for k in data["seq"].keys():
        search[k]=[x*ws for x in search[k]]        
        
    #geo filtering
    if(geo!="none"):
        if(geo!="intergenic"):
            for k in data["seq"].keys():
                sk=search[k]
                
                if(geo!="RNA_gene"):
                    tt=[geo]
                else:
                    tt=["ncRNA_gene", "tRNA_gene", "snRNA_gene", "snoRNA_gene", "rRNA_gene"]
                if(len(sk)>0):
                    #p=[x*ws for x in sk]                    
                    annot=annotationsLocal(positions=sk, window=ws*len(pattern), gff=data["gff"][k], types=tt, onlyIDs="False", intersect=intersect)
                    #p=[x/ws for x in annot.keys()];
                    p=[(int)(x) for x in annot.keys()];
                    search[k]=p
                else:
                    search[k]=[]
        else:   #intergenic regions, by now without any annotations
            for k in data["seq"].keys():
                sk=search[k]
                ws=data["windowSize"]
                
                if(len(sk)>0):
                    #p=[x*ws for x in sk]
                    
                    annot=annotationsLocal(positions=sk, window=ws*len(pattern), gff=data["gff"][k], types=["gene", "pseudogene", "ncRNA_gene", "tRNA_gene", "snoRNA_gene", "snRNA_gene"], onlyIDs="False", intersect=intersect)
                    p1=set(p)-set(annot.keys())
                    #p=[x/ws for x in p1];
                    p=[(int)(x) for x in p1];
                    search[k]=list(p)
                else:
                    search[k]=[]
       
    #for JSON serialization   
    for k in search.keys():
        search[k]=(str)(search[k])
            
    data["search"]={'points':search, 'sizePattern':len(pattern)}
    return jsonify(points=search, sizePattern=len(pattern))



#%%
@app.route("/getPartSeq")
def getPartSeq(start=0, end=0, track="None"):
    global data
    import numpy as np
    start=int(request.args.get("start"))
    end=int(request.args.get("end"))
    track=str(request.args.get("track"))
    seq=data["seq"][track]
    part=np.array(seq[start:end],dtype=float)
    return jsonify(partSeq=list(part))



#%%
@app.route("/annotations")
def annotations(positions=[], window=1000, types=["any"], track="None", onlyIDs="False", align="left", intersect="soft"):
    window=int(request.args.get("window"))
    pos=eval(request.args.get("positions"))
    types=eval(request.args.get("types"))
    track=str(request.args.get("track"))
    onlyIDs=str(request.args.get("onlyIDs"))
    align=str(request.args.get("align"))
    intersect=str(request.args.get("intersect"))

    print("INTERSECT ANNOT is",intersect)
    
    res=annotationsLocal(positions=pos, gff=data["gff"][track], window=window, types=types, onlyIDs=onlyIDs, align=align, intersect=intersect)
    return jsonify(response=res)

def annotationsLocal(positions, gff, window=1000, types=["any"], track="None", onlyIDs="False", align="left", intersect="soft"):
    import time
    global data
    t0=time.clock()
    res=ann.annotate(positions, gff, types, window, align, intersect)
    print("Annotations take",(time.clock()-t0),"s")
    if(onlyIDs=="True"):
        ids=[]
        for k in res.keys():
            kl=res[k]
            for x in kl:
                ids.append(x["id"])
        return list(set(ids))
    return res


#%%
@app.route("/nucleotides")
def nucleotides(start=0, end=10, track="None"):
    global data

    start=int(request.args.get("start"))
    end=int(request.args.get("end"))
    track=str(request.args.get("track"))
    print(data["fasta"].keys())
    return jsonify(response=data["fasta"][track][start:end])
 
#%%
'''
Given a number of genomic positions on a given wig-fasta track, 
it returns the nucleotide profile for such positions
Size refers to the sequence length to take for computing the profiles, 
starting at positions.

Originally, this method aligned profiles but that was very costly computationally
(commented code). Now it only retrieves the sequences and computes the best
k-motif (currently also limited to 6-mers)
'''
@app.route("/nucProfile")
def nucProfile(positions=[], size=10, track="None"):
    global data
    import time
    t00=time.clock()
    pos=eval(request.args.get("positions"))
    size=int(request.args.get("size"))
    track=str(request.args.get("track"))
    #align=str(request.args.get("align"))
    
    #basic operations (prof and c might not be very useful)
    seqs={}
    for p in pos:
        seqs[p]=data["fasta"][track][p:p+size].upper()
    c=ms.consensus(seqs.values())
    prof=ms.profile(seqs.values())

    #gibbs sampling for motif search
    if(len(seqs)>1):
        t0=time.clock()
        for x in range(6,7,1): #for different k
            for i in range(1):
                bm=ms.gibbsSampler(seqs.values(),x,1000)
                print(i,")",x, bm["score"], bm["score"]/len(seqs))
        print("Whole Gibbs took",(time.clock()-t0))    
    else:
        bm={}
    mot={}
    motloc={}
    for i in range(len(bm["motifs"])):
               mot[seqs.keys()[i]]=bm["motifs"][i]
               motloc[seqs.keys()[i]]=seqs.values()[i].find(bm["motifs"][i])
#           
#    #Alignment is superslow  we only do it if #seqx<50 with kalign to keep times <=1s
#    if(len(pos)<=50):
#        align=helpers.align(seqs, method="kalign_msa")
#        ac=ms.consensus(align.values())
#        aprof=ms.profile(align.values())
#        print("nucProfile with alignment took",(time.clock()-t00))
#        return jsonify(seqs=seqs, consensus=c, alignment=align, aconsensus=ac, aprofile=aprof)
#    print("No alignment performed: too many sequences")
    print("nucProfile took",(time.clock()-t00))
    return jsonify(seqs=seqs, consensus=c, profile=prof, motifs=mot, locations=motloc, motifConsensus=bm["consensus"], motifProfile=bm["profile"])
 
#%%
@app.route("/annotationsGOA")
#annotations is a list of gene IDs? TODO: reshaping
#types by now is any. Maybe in the future we divide in BP,MF,CC
def annotationsGOA(genes=[], types=["any"]):
    global data
    
    a=eval(request.args.get("genes"))
    print('annotation size: ',len(a))
    agon=ann.annotateGOnames(a, data["goa"], data["go"])
    return jsonify(response=agon)
    
#%%
@app.route("/enrichmentGO")
#annotations is the result of calling annotations
#correction is for multiple hipothesis and can be 'none', 'bonferroni', 'fdr' or 'fwer'
#alpha is the threshold, applied as is with 'none' correction or with the specified one
def enrichmentGO(annotations={}, correction="none", alpha=0.01):
    global data
    print("Enrichment GO")
    gis=set(eval(request.args.get("annotations")))
    alpha=float(request.args.get("alpha"))
    correction=request.args.get("correction")
    
    #print gis
#    gis=set()
#    for x in annotations.keys(): #TODO here
#        for y in annotations[x]:
#            gis.add(y["id"])
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

#@app.after_request
#def serialize_passport(response):
#    if hasattr(g, "passport"):
#        session["passport_id"] = g.passport.id
#    return response
    
@app.route("/test")
def test():
    d={"a":[1,2,3], "b":[4,5,6]}
    return jsonify(d)

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
    #app.run(debug=True)
    app.run(debug=True, host='0.0.0.0', port=2750)
    
    
    