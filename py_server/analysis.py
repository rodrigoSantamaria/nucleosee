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
    
global session
global data
session={}
data={}

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
           
           
@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    global user
    import time
    
    t0=time.clock()
    print("Uploading file...")
    if request.method == 'POST':
        print("requesting...")
        file = request.files['file']#this is the data transmission
        print("file requested")
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            root=os.path.join(app.config['UPLOAD_FOLDER'])
            #We decided not to use any 'user' infraestructure
            #if (user in os.listdir(root))==False:                 #create directory for this user
            #    os.mkdir(os.path.join(root,user))
            #if filename in os.listdir(os.path.join(root,user)):
            if filename in os.listdir(root):
                print(filename,'already uploaded')    #TODO: by now we avoid resubmissions (for tests)
            else:
                print('uploading...')
                file.save(os.path.join(root, filename))
                #file.save(os.path.join(root, user, filename))
                print("file saved in ",(time.clock()-t0),"s")
            #return redirect(url_for('uploaded_file', filename=filename))
            return jsonify(path=os.path.join(app.config['UPLOAD_FOLDER'], filename))
        else:
            print("file not found", file.filename, allowed_file(file.filename))
            return json_response(status=415, responseText="File does not exist or unsported file type (.wig or .bw allowed)")
#%%
from flask import send_from_directory

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)
      
#%%
@app.route('/testUpload', methods=['GET'])
def testUpload():
    import re
    filename=""+request.args.get("filename") 
    print("Testing upload of ",filename)
    cpath=os.path.join(app.config['UPLOAD_FOLDER'],user)
    if os.path.exists(cpath) and filename in os.listdir(cpath):
        #.wig exists
        filename=re.sub(r"\..*$", "", filename)
        for f in os.listdir(cpath):
            try:
                if(f.startswith(filename) and f.endswith(".pic")):
                    #.pic exists
                    f=re.sub("^"+filename,"",f)
                    org=re.sub(".*org", "", re.sub("\\.pic$","",f))
                    f=re.sub("org"+org, "",f)
                    org=re.sub("_", " ", org)
                    clip=float(re.sub("ws.*$","",re.sub("c","",f)))
                    f=re.sub("c"+str(clip), "", f)
                    ws=int(re.sub("nb.*$","",re.sub("ws","",f)))
                    f=re.sub("ws"+str(ws), "", f)
                    nb=int(re.sub("i.*$","",re.sub("nb","",f)))
                    f=re.sub("nb"+str(nb), "", f)
                    interpol=re.sub(".pic$", "", re.sub("^i","",f))
                    
                    #print(clip," ",ws," ",nb, " ", interpol, " ", org)
                    return jsonify(response="file exists", clip=clip, ws=ws, nb=nb, interpol=interpol, org=org)
            except:
                continue
        return jsonify(response='outdated version') #no .pic (preprocess required)
    else:
        return jsonify(response='not found')#no .wig (upload required)
#testUpload("dwtMini2.wig")
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
@app.route("/loadData")
def loadData(dataName="Test", track="None", clear="false"):
    import time
    import os
    import annotations as ann
    import pickle #i tried cPickle but is way slower!
    global session
    global user
    
    dataName=request.args.get("dataName")
    t00=time.clock()
    clear=request.args.get("clear")
    if (clear=="true"):
        session[user]={}
    
    
    #0) Get info from the track file
    #basePath=os.path.join(app.config['UPLOAD_FOLDER'],user)
    basePath=app.config['UPLOAD_FOLDER']
    f=open(os.path.join(basePath,"tracks.txt"))
    for l in f.readlines():
        if(l.split("\t")[0].strip()==dataName):
            break

    chars=l.split("\t") 
    picklePaths=chars[1].strip().split(",")
    if(len(picklePaths)>1):
        picklePath=chars[0].strip()+".pic"
    else:
        picklePath=picklePaths[0]
   
    organism=chars[4].strip()
    
    #----Load the corresponding pickle
    data={}
    
    t0=time.clock()
    print("---------------- FILE: "+picklePath+" ------------------")
    f=open(os.path.join(basePath,picklePath))
    pdata=pickle.load(f)
    print("KEYS", pdata.keys())
    
    filenames=filter(lambda x:x.endswith(".wig") or x.endswith(".bw"),pdata.keys())
    if(len(filenames)>0):
        filename=filenames[0]
    else: #batch data
        filename="processed"
 
    if track=="None":
        track=sorted(pdata[filename]["maximum"].keys())[0]
    if(len(pdata[filename]["gff"].keys())!=len(pdata[filename]["seq"].keys())):
            print("One or more chromosome tracks do not match with GFF names:")      
    
   
    print("load data takes",(time.clock()-t0))
    
    #--------------- COMPILE BATCH
    data["batch"]={}
    data["batch"]["processed"]=pdata[filename]
    if(len(filenames)==1): #  (single file)
        data["batch"]["min"]=pdata[filename]["seq"]
        data["batch"]["max"]=pdata[filename]["seq"]
        data["go"]=pdata["go"]
        data["goa"]=pdata["goa"]
        data["windowSize"]=pdata["windowSize"]
    else:               # multiple files
        data["batch"]["min"]=pdata["max"]
        data["batch"]["max"]=pdata["min"]
        
        #here we require to get the go/goa data (can be fixed)
        t0=time.clock()
        dataGO=ann.go()
        dataGOA=[]
    
        try:
            dataGOA=ann.goa(organism)
            print("GO loaded")
        except:
            print("organism",organism,"'s annotations missing or failing")
        print("done! ... GO annotations takes",(time.clock()-t0))
        
        data["go"]=dataGO
        data["goa"]=dataGOA       
        data["windowSize"]=pdata["processed"]["windowSize"]
    
  
    session[user][dataName]=data#trying to load several datasets at once
    print('file preprocess takes ',(time.clock()-t00),"s")
    print("USER KEYS:", session[user].keys())
    
    dbp=data["batch"]["processed"]
    
    return jsonify(seq=dbp["res"][track], 
                   fullLength=dbp["fullLength"],
                   maximum=(float)(dbp["maximum"][track]),
                   minimum=(float)(dbp["minimum"][track]), 
                    mean=(float)(dbp["mean"][track]), 
                stdev=(float)(dbp["stdev"][track]), 
                dseq=dbp["dseq"][track], 
                chromosomes=sorted(dbp["maximum"].keys()),
                bins=dbp["bins"][track],
                windowSize=data["windowSize"])

#ret=loadData("Test")
#%%
@app.route("/listData")
def listData():
    try:
        f=open(os.path.join(app.config['UPLOAD_FOLDER'],"tracks.txt"))
    except:
        return jsonify(response="error", msg="no tracks.txt at {} in {}".format(app.config['UPLOAD_FOLDER'],os.getcwd()))
    tracks=[l.split("\t")[0] for l in f.readlines()]
    return jsonify(response=tracks)

#%%
"""
Preprocess as wig or bw file
filenames  array with file names (to be found in the user folder under UPLOAD_FOLDER)
dataName   name of the resulting preprocessed data (def "None")
windowSize size of the discretization window (def 100)
numBins    number of discretization bins (def 5). 
            E.g. if windowSize=100 and numBins=5, the method takes 100 nucleotides
            computes their average level and assign them a discrete level out of
            5 possible levels (which in the default mode implies 5 percentile ranges)
maxSize    maximum length of the returning array (to avoid bandwidth overload). Default 100K
            It will sample the data to fit this value. E.g. if maxSize=100K and
            total size is 200K it will return an array with the average level of
            every 2 original values.
stdev      outlier clipping. It will clip levels above/below this number of 
            standard deviations from the mean. Default 3
track      in the case of several tracks (usually chromosomes), which one to return. 
            Default 'None' returns the first track on alphabetic order
organism   species name for the organism (standard names such as 'Saccharomyces cerevisiae')
            It will be used for annotations' (GFF, GOA, FASTA) preloading
interpolation In the case of bigWig files with missing (variableStep) values,
            this value indicates how to impute them. Defalult "mean"
returns    a JSON object with the following fields:
               result   valores normalizados (sin discretizar) sampleados hasta maxSize
               fullLength  longitud total de los datos iniciales/normaliados
               maximum,minimum,mean,sdev   medidas estadísticas de los datos normalizados
               dseq        datos discretizados por la función discretize()
"""
@app.route("/preprocess")
def batchPreprocess(filenames=[], dataName="None", windowSize=100, numBins=5, maxSize=100000, stdev=3, track="None", organism="Saccharomyces cerevisiae", interpolation="mean"):
    global data
    global session
#
#    
    import time
    import re
    import os, pickle        
    import numpy as np
    import helpers

    t00=time.clock()
#    
#    #0) getting parameter
    print('----------------------- PREPROCESS...')
    t0=time.clock()
    basePath=app.config['UPLOAD_FOLDER']
    filenames=eval(request.args.get("filenames"))
    track=request.args.get("track")
    organism=request.args.get("organism")
    windowSize=int(request.args.get("windowSize"))
    numBins=int(request.args.get("numBins"))
    maxSize=int(request.args.get("maxSize"))
    interpolation=request.args.get("interpolation")
    stdev=float(request.args.get("stdev"))
    dataName=request.args.get("dataName")
    clear=request.args.get("clear")

    data={}
    print("Files are: ",filenames)
    
    #TODO: check that all batch files have same tracks and sizes
    #3) annotations (same for all batch files? - hard to do with pickles)
    t0=time.clock()
    dataGO=ann.go()
    dataGOA=[]

    try:
        dataGOA=ann.goa(organism)
        print("GO loaded")
    except:
        print("organism",organism,"'s annotations missing or failing")
    print("done! ... GO annotations takes",(time.clock()-t0))
    
    data["go"]=dataGO
    data["goa"]=dataGOA
    data["windowSize"]=windowSize

    #basePath="/Users/rodri/WebstormProjects/seqview/py_server/genomes/jpiriz/"
    print("organism:",organism)
        
    pickleFiles=[]
    #------------ Preprocess each file
    for filename in filenames: #for each file
        pickleFile=re.sub(r"\..*$", "", filename)+"c"+str(stdev)+"ws"+str(windowSize)+"nb"+str(numBins)+"i"+interpolation+"org"+organism.replace(" ", "_")+".pic"
        picklePath=os.path.join(basePath,pickleFile)
        path=os.path.join(basePath,filename)
        pickleFiles.append(pickleFile)

        savePickle=False
        
        print("---------------- FILE: "+filename+" ------------------")
        genome={}
        #1A) Preprocessed data (.pic data) exist
        if os.path.isfile(picklePath):
            f=open(picklePath)
            data[filename]=pickle.load(f)
            print("pickle loaded!", data[filename].keys())
            print(data.keys())
            
            talg3=data[filename]["gff"]["chromosome3"]
            print("CDS MATCH: ", talg3[talg3["s"]==115781])
            
            if track=="None":
                track=sorted(data[filename]["maximum"].keys())[0]
                    
        #1B) Preprocessing must be done     
        else: #no previous preprocessing -> do it now
            if(path.endswith("bw")):
                tbw=time.clock()
                data[filename]=helpers.processBigWig(path=path, track=track, windowSize=windowSize, numBins=numBins, percentile=True, maxSize=maxSize, stdev=stdev, organism=organism, picklePath=picklePath)
                print("process bw done in ",(time.clock()-tbw))
                genome=data[filename]["seq"]
    
            else:
                genome=helpers.readWig(path, method=interpolation)
                data[filename]=helpers.processWig(genome, stdev=stdev, windowSize=windowSize, numBins=numBins, maxSize=maxSize, percentile=True, organism=organism, track=track)
    
            if track=="None":
                track=sorted(genome.keys())[0]
            savePickle=True
    
        
        print("load data takes",(time.clock()-t0))
       
        if(len(data[filename]["gff"].keys())!=len(data[filename]["seq"].keys())):
            print("One or more chromosome tracks do not match with GFF names:")
            print(data[filename]["gff"].keys())
            print(data[filename]["seq"].keys())
            
        if(savePickle):
            tpickle=time.clock()
            print("saving pickle in path:", picklePath)
            print("with keys",data.keys())
            print("with keys",data[filename].keys())
            f=open(picklePath, 'w')
            pickle.dump(data,f)
            f.close()
     
    print("Compiling batch")  
    #--------------- COMPILE BATCH
    data["batch"]={"min":{}, "max":{}}
    meanBatch={}
    if(len(filenames)>1):
        if(dataName+".pic" in os.listdir(basePath))==False:
            print("Computing and saving batch data")
            tpickle=time.clock()
            for k in data[filenames[0]]["seq"].keys():
                data["batch"]["min"][k]=np.minimum(data[filenames[0]]["seq"][k], data[filenames[1]]["seq"][k])
                data["batch"]["max"][k]=np.maximum(data[filenames[0]]["seq"][k], data[filenames[1]]["seq"][k])
                temp=np.ndarray(shape=(len(filenames), len(data[filenames[0]]["seq"][k])))
                temp[0]=data[filenames[0]]["seq"][k]
                temp[1]=data[filenames[1]]["seq"][k]
                for i in range(2,len(filenames)):
                    data["batch"]["min"][k]=np.minimum(data["batch"]["min"][k], data[filenames[i]]["seq"][k])
                    data["batch"]["max"][k]=np.maximum(data["batch"]["max"][k], data[filenames[i]]["seq"][k])
                    temp[i]=data[filenames[i]]["seq"][k]
                meanBatch[k]=np.mean(temp,axis=0) 
            data["batch"]["processed"]=helpers.processWig(meanBatch, stdev, windowSize, numBins, maxSize, True, organism, k)
            f=open(os.path.join(basePath,dataName+".pic"), 'w')
            pickle.dump(data["batch"],f)
            f.close()
            print("serialize data takes",(time.clock()-tpickle))
        else:
            print("Loading batch data")
            f=open(os.path.join(basePath,dataName+".pic"), 'r')
            data["batch"]=pickle.load(f)
            
    else:
        data["batch"]["processed"]=data[filenames[0]]
        data["batch"]["min"]=data[filenames[0]]["seq"]
        data["batch"]["max"]=data[filenames[0]]["seq"]
    
    #add entry on listing file 
    if(savePickle):
        f=open(os.path.join(basePath,"tracks.txt"), 'a')
        f.write(dataName+"\t"+",".join(pickleFiles)+"\t"+str(windowSize)+"\t"+str(numBins)+"\t"+organism+"\n")
        f.close()
    
    session[user]=data
    print('FILE PREPROCESS TAKES ',(time.clock()-t00),"s")
    print("Data KEYS: ",data.keys())
    
    return jsonify(seq=data["batch"]["processed"]["res"][track], 
                   fullLength=data["batch"]["processed"]["fullLength"],#len(genome[track]), 
                   maximum=(float)(data["batch"]["processed"]["maximum"][track]),
                   minimum=(float)(data["batch"]["processed"]["minimum"][track]), 
                    mean=(float)(data["batch"]["processed"]["mean"][track]), 
                stdev=(float)(data["batch"]["processed"]["stdev"][track]), 
                dseq=data["batch"]["processed"]["dseq"][track], 
                chromosomes=sorted(data["batch"]["processed"]["maximum"].keys()),
                bins=data["batch"]["processed"]["bins"][track])
#    return data
   
#batch=["MN-Damien-WT-1_S1_wlt_mean.wig", "MN-Damien-WT-2_S2_wlt_mean.wig"]                
#data=batchPreprocess(filenames=batch, windowSize=100, numBins=3, maxSize=100000, stdev=3, track="None", organism="Schizosaccharomyces pombe", interpolation="mean")
     
#%%
#data["batch"]={"min":{}, "max":{}, "mean":{}}
#for k in data[batch[0]]["seq"].keys():
#    data["batch"]["min"][k]=np.minimum(data[batch[0]]["seq"][k], data[batch[1]]["seq"][k])
#    data["batch"]["max"][k]=np.maximum(data[batch[0]]["seq"][k], data[batch[1]]["seq"][k])
#    temp=np.ndarray(shape=(len(batch), len(data[batch[0]]["seq"][k])))
#    temp[0]=data[batch[0]]["seq"][k]
#    temp[1]=data[batch[1]]["seq"][k]
#    for i in range(2,len(batch)):
#        data["batch"]["min"][k]=np.minimum(data["batch"]["min"][k], data[batch[i]]["seq"][k])
#        data["batch"]["max"][k]=np.maximum(data["batch"]["max"][k], data[batch[i]]["seq"][k])
#        temp[i]=data[batch[i]]["seq"][k]
#    data["batch"]["mean"][k]=np.mean(temp,axis=0)  
#data["batch"]["processed"]=helpers.processWig(data["batch"]["mean"], stdev=3, windowSize=100, numBins=3, maxSize=100000, percentile=True, organism="Schizosaccharomyces pombe", track=k)
    
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
#%%#
def preprocess(filename="dwtMini2.wig", windowSize=100, numBins=5, maxSize=100000, stdev=3, track="None", recharge="False", organism="Saccharomyces cerevisiae", interpolation="mean"):
    global data
    global session

    
    t00=time.clock()
    
    #0) getting parameter
    print('Parameters...')
    t0=time.clock()
#    basePath=os.path.join(app.config['UPLOAD_FOLDER'],user)
    basePath=app.config['UPLOAD_FOLDER']
    filename=str(request.args.get("filename"))
    path=os.path.join(basePath,filename)
    track=request.args.get("track")
    forceReload=request.args.get("recharge")
    organism=request.args.get("organism")
    windowSize=int(request.args.get("windowSize"))
    numBins=int(request.args.get("numBins"))
    maxSize=int(request.args.get("maxSize"))
    interpolation=request.args.get("interpolation")
    stdev=float(request.args.get("stdev"))
    f=re.sub(r"\..*$", "", filename)+"c"+str(stdev)+"ws"+str(windowSize)+"nb"+str(numBins)+"i"+interpolation+"org"+organism.replace(" ", "_")+".pic"
    picklePath=os.path.join(basePath,f)
    savePickle=False
    
    print("File: "+f)
    print("organism:",organism)
    print("reload:", forceReload)
    
    genome={}
    
    print("---> picklePath", picklePath, os.path.isfile(picklePath))
    
    #1A) Preprocessed data (.pic data) exist
    if os.path.isfile(picklePath) and forceReload=="false":
        print("Pickle exists!!!, recharge=", forceReload)
        f=open(picklePath)
        data=pickle.load(f)
        print("pickle loaded!", data["maximum"].keys())
        
        if track=="None":
            track=sorted(data["maximum"].keys())[0]
            
        print ("track is ", track)

    #1B) Preprocessing must be done     
    else: #no previous preprocessing -> do it now
        if(path.endswith("bw")):
            tbw=time.clock()
            data=helpers.processBigWig(path=path, track=track, windowSize=windowSize, numBins=numBins, percentile=True, maxSize=maxSize, stdev=stdev, organism=organism, picklePath=picklePath)
            print("process bw done in ",(time.clock()-tbw))
            genome=data["seq"]

        else:
            genome=helpers.readWig(path, method=interpolation)
            data=helpers.processWig(genome, stdev=stdev, windowSize=windowSize, numBins=numBins, maxSize=maxSize, percentile=True, organism=organism, track=track)

        if track=="None":
            track=sorted(genome.keys())[0]
        savePickle=True

    
    print("load data takes",(time.clock()-t0))
   
    if(len(data["gff"].keys())!=len(data["seq"].keys())):
        print("One or more chromosome tracks do not match with GFF names:")
    print("\tData names:",data["seq"].keys())
    print("\tGFF names:", data["gff"].keys())

    #3) annotations
    t0=time.clock()
    dataGO=ann.go()
    dataGOA=[]

    try:
        dataGOA=ann.goa(organism)
        print("GO loaded")
    except:
        print("organism",organism,"'s annotations missing or failing")
    print("done! ... GO annotations takes",(time.clock()-t0))
    
    data["go"]=dataGO
    data["goa"]=dataGOA
    data["windowSize"]=windowSize
       
    print("Setting session data")
    session[user]=data
    print("Save pickle?", savePickle)
    
    if(savePickle):
        tpickle=time.clock()
        print("saving pickle in path:", picklePath)
        datap=data

        f=open(picklePath, 'w')
        pickle.dump(datap,f)
        print("serialize data takes",(time.clock()-tpickle))
        #NOTE: should remove the .wig here to avoid double memory space?

    print('whole preprocess takes ',(time.clock()-t00),"s")
    
    print('returning track',track)
    print('bins are',data["bins"][track])
    #data["bins"][track]=[abs(x) for x in data["bins"][track]]
    
    return jsonify(seq=data["res"][track], 
                   fullLength=data["fullLength"],#len(genome[track]), 
                   maximum=(float)(data["maximum"][track]),
                   minimum=(float)(data["minimum"][track]), 
                    mean=(float)(data["mean"][track]), 
                stdev=(float)(data["stdev"][track]), 
                dseq=data["dseq"][track], 
                chromosomes=sorted(data["maximum"].keys()),
                bins=data["bins"][track])

#%%preprocess(filename="/Users/rodri/Documents/investigacion/IBFG/nucleosomas/dwtMini2.wig")

#%%
@app.route("/getTrack")
def getTrack(track="None", dataName="None"):
    global data
    
    track=request.args.get("track")
    dataName=request.args.get("dataName")
    print("returning chromosome",track," for ",dataName)
    if track=="None":
        track=data["res"].keys()[0]

    if( ("search" in data.keys()) == False):
        data["search"]={"points":{}, "sizePattern":0}
    if( ("ego" in data.keys()) == False):
        data["ego"]={}
    
#    return jsonify(seq=data["batch"]["processed"]["res"][track], fullLength=len(data["batch"]["processed"]["seq"][track]), maximum=(float)(data["batch"]["processed"]["maximum"][track]), minimum=(float)(data["batch"]["processed"]["minimum"][track]), mean=(float)(data["batch"]["processed"]["mean"][track]), stdev=(float)(data["batch"]["processed"]["stdev"][track]), dseq=data["dseq"][track], bins=data["batch"]["processed"]["bins"][track], chromosomes=sorted(data["batch"]["processed"]["res"].keys()), search=data["search"], ego=data["ego"])
    if(dataName=="None"):
        dbp=data["batch"]["processed"]
    else:
        dbp=data[dataName]["batch"]["processed"]
        
    return jsonify(seq=dbp["res"][track], fullLength=dbp["fullLength"],
            maximum=(float)(dbp["maximum"][track]), minimum=(float)(dbp["minimum"][track]), 
            mean=(float)(dbp["mean"][track]), stdev=(float)(dbp["stdev"][track]),
             dseq=dbp["dseq"][track], bins=dbp["bins"][track], chromosomes=sorted(dbp["res"].keys()),
             search=data["search"], ego=data["ego"])

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
def search(pattern="", d=0, geo="none", intersect="soft", softMutations="false", 
           dataName1="None", dataName2="None", pattern2="None", join="not"):
    global data
    data["search"]={}
    data["gis"]=set()
    data["annotations"]={}
    
    print("searching..:")
    print(data.keys())
    
    d=int(request.args.get("d"))
    pattern=str(request.args.get("pattern"))
    geo=str(request.args.get("geo"))
    intersect=str(request.args.get("intersect"))
    softMutations=str(request.args.get("softMutations"))
    dataName1=str(request.args.get("dataName1"))
    dataName2=str(request.args.get("dataName2"))
    join=str(request.args.get("join"))
    pattern2=str(request.args.get("pattern2"))

    print("SEARCHING IN DATANAME1", dataName1)   
    if dataName1=="None":
        ret=searchLocal(data, pattern, d, geo, intersect, softMutations)
    else:
        ret=searchLocal(data[dataName1], pattern, d, geo, intersect, softMutations)
    search1=ret["points"]
    #sizePattern=ret["sizePattern"]
    
    if(dataName2=="None" or join=="None"):
        search=search1
    else:
        if(pattern2=="None"):
            ret=searchLocal(data[dataName2], pattern, d, geo, intersect, softMutations)
        else:
            ret=searchLocal(data[dataName2], pattern2, d, geo, intersect, softMutations)
            
        search2=ret["points"]
        search={}
        for k in search2.keys():
            s1=set(search1[k])
            s2=set(search2[k])
            if(join=="not"):
                s12=s1 - s2
            if(join=="and"):
                s12=s1 & s2
            if(join=="or"):
                s12=s1 | s2
                if(len(data["search"].keys())==0):
                    data["search"][dataName1]={"points":{}, "sizePattern":0}
                    data["search"][dataName2]={"points":{}, "sizePattern":0}
                data["search"][dataName1]["points"][k]=list(s1)
                data["search"][dataName2]["points"][k]=list(s2)
            search[k]=list(s12)
            print("SEARCH SIZES: ", len(s1), " ", len(s2), " ", len(s12))    
    print("DATA NAMES:" ,dataName1, dataName2)
    
    if(join=="None" or dataName2=="None" or join=="and" or join=="not"):
        data["search"][dataName1]={"points":search,"sizePattern":ret["sizePattern"]}#Asume same len pattern on multi-searches
    if(join=="and" and dataName2!="None"):
        data["search"][dataName2]={"points":search,"sizePattern":ret["sizePattern"]}#Asume same len pattern on multi-searches
    if(join=="or" and dataName2!="None"):
        data["search"][dataName1]["sizePattern"]=ret["sizePattern"]
        data["search"][dataName2]["sizePattern"]=ret["sizePattern"]
    #data["search"]={"points":search,"sizePattern":ret["sizePattern"]}#Asume same len pattern on multi-searches
    
    search=data["search"]
    #for json
    for k in search.keys():
        #search[k]=(str)(search[k])
        for j in search[k]["points"]:
            search[k]["points"][j]=(str)(search[k]["points"][j])
        
    if("response" in ret.keys()):
        return jsonify(response=ret["response"], msg=ret["msg"])
    else:
        #return jsonify(points=search, sizePattern=sizePattern)
        return jsonify(search)
        
#%%    
def searchLocal(data, pattern="", d=0, geo="none", intersect="soft", softMutations="false"):
    import time
    
    t00=time.clock()
    ws=data["windowSize"]
                
    #CASE 1) Numerical range pattern
    if(string.find(pattern, "go:")==-1):#go terms might contain '-'
        interval=helpers.convertRange(pattern)
        if(interval!=-1):
            print("NUMERICAL RANGE")
            points={}
            for k in data["batch"]["processed"]["seq"].keys():
                points[k]=(str)([(int)(interval["start"])] if (interval["start"]+interval["length"])<len(data["batch"]["processed"]["seq"][k]) else [])
            data["search"]={'points':points, 'sizePattern':((int)(interval["length"]/ws))}
            return {"points":points, "sizePattern":((int)(interval["length"]/ws))}
    
    #CASE 2) gene/go name pattern
    t=data["batch"]["processed"]["bwt"][data["batch"]["processed"]["seq"].keys()[0]]
    
    pattern=pattern.lower().strip()
    patternLetters=[chr(x) for x in range(ord('a'), ord('a')+len(data["batch"]["processed"]["bins"][data["batch"]["processed"]["bins"].keys()[0]])-1)]

    patternSymbols=patternLetters
    patternSymbols.append('+')
    patternSymbols.append('*')
    for x in range(0,10):
        patternSymbols.append((str)(x)) 
    
    print(patternSymbols, "\t",set(pattern))
    if(False in [x in patternSymbols for x in set(pattern)]):
        print("GENE OR TERM")      
        points={}
        sizes={}
        for k in data["batch"]["processed"]["seq"].keys():
            if(string.find(pattern,"go:")==-1):
                print("GENE")
                oc=ann.searchGene(pattern, data["batch"]["processed"]["gff"][k], ["gene", "tRNA_gene"])
            else:
                print("TERM")
                oc=ann.searchGO(pattern.replace("go:", "").strip(), data["goa"], data["go"], data["batch"]["processed"]["gff"][k])
            points[k]=(str)([(int)(y["start"]) for y in oc])
            sizes[k]=(str)([(int)((y["end"]-y["start"])/ws) for y in oc])
        data["search"]={'points':points, 'sizePattern':sizes}
        #return jsonify(points=points, sizePattern=sizes);
        return {"points":points, "sizePattern":sizes};
        
    #CASE 3) bin pattern    
    pattern=helpers.convertString(pattern)
    print("BWT ON ",pattern)
    search={}
    for k in data["batch"]["processed"]["seq"].keys():
         if(False in [x in patternLetters for x in set(pattern)]):
            search[k]=[] #instead of returning an error
         else:
            t0=time.clock()
            t=data["batch"]["processed"]["bwt"][k]
            search[k]=(ss.bwMatchingV8("".join(data["batch"]["processed"]["dseq"][k]), pattern, t["bwt"], t["firstOccurrence"],t["suffixArray"],t["checkpoints"],1000, d))
            print(len("".join(data["batch"]["processed"]["dseq"][k])))
            print("Search ",k,"takes",(time.clock()-t0), "and finds",len(search[k]), "occurences")
            if(len(search[k])>10000):
                return {"response":"error", "msg":"Too many occurrences, please narrow your search", "points":{}, "sizePattern":len(pattern)}
    print("Search finished in ",(time.clock()-t00))
    
    
    if(softMutations=="true"):
        ti=time.clock()
        for k in data["batch"]["processed"]["seq"].keys():
            search[k]=helpers.filterHard(search[k], data["batch"]["processed"]["dseq"][k], pattern)
        print("Soft mutation filtering takes ", (time.clock()-ti))
    
    for k in data["batch"]["processed"]["seq"].keys():
        search[k]=[x*ws for x in search[k]]        
        
    #geo filtering
    if(geo!="none"):
        if(geo!="intergenic"):
            for k in data["batch"]["processed"]["seq"].keys():
                sk=search[k]
                
                if(geo!="RNA_gene"):
                    tt=[geo]
                else:
                    tt=["ncRNA_gene", "tRNA_gene", "snRNA_gene", "snoRNA_gene", "rRNA_gene"]
                if(len(sk)>0 and (k in data["batch"]["processed"]["gff"].keys())):
                    annot=annotationsLocal(positions=sk, window=ws*len(pattern), gff=data["batch"]["processed"]["gff"][k], types=tt, onlyIDs="False", intersect=intersect)
                    p=[(int)(x) for x in annot.keys()];
                    search[k]=p
                else:
                    search[k]=[]
        else:   #intergenic regions, by now without any annotations
            for k in data["batch"]["processed"]["seq"].keys():
                sk=search[k]
                ws=data["windowSize"]
                
                if(len(sk)>0 and (k in data["batch"]["processed"]["gff"].keys())):
                    if(intersect=="soft"):
                        intersect2="hard"
                    else:
                        intersect2="soft"
                    annot=annotationsLocal(positions=sk, window=ws*len(pattern), gff=data["batch"]["processed"]["gff"][k], types=["gene", "pseudogene", "ncRNA_gene", "tRNA_gene", "snoRNA_gene", "snRNA_gene"], onlyIDs="False", intersect=intersect2)
                    p1=set(sk)-set(annot.keys())
                    p=[(int)(x) for x in p1];
                    search[k]=list(p)
                else:
                    search[k]=[]
       
    return {"points":search, "sizePattern":len(pattern)}



#%%
@app.route("/getPartSeq")
def getPartSeq(start=0, end=0, track="None", dataName="None"):
    global data
    import numpy as np
    start=int(request.args.get("start"))
    end=int(request.args.get("end"))
    track=str(request.args.get("track"))
    dataName=str(request.args.get("dataName"))
    
    if(dataName=="None"):
        seq=data["batch"]["processed"]["seq"][track]
        part=np.array(seq[start:end],dtype=float)
        partMax=np.array(data["batch"]["max"][track][start:end],dtype=float)
        partMin=np.array(data["batch"]["min"][track][start:end],dtype=float)
    else:
        seq=data[dataName]["batch"]["processed"]["seq"][track]
        part=np.array(seq[start:end],dtype=float)
        partMax=np.array(data[dataName]["batch"]["max"][track][start:end],dtype=float)
        partMin=np.array(data[dataName]["batch"]["min"][track][start:end],dtype=float)
     
    return jsonify(partSeq=list(part), minimum=list(partMin), maximum=list(partMax))

#%%
@app.route("/getDSeq")
def getDSeq(start=0, end=0, track="None", dataName="None"):
    global data
   
    start=int(request.args.get("start"))
    end=int(request.args.get("end"))
    track=str(request.args.get("track"))
    dataName=str(request.args.get("dataName"))
    
    if(dataName=="None"):
        seq=data["batch"]["processed"]["dseq"][track]
        ws=data["batch"]["processed"]["windowSize"]
        part=seq[int(start/ws):int(end/ws)]
    else:
        seq=data[dataName]["batch"]["processed"]["dseq"][track]
        ws=data[dataName]["batch"]["processed"]["windowSize"]
        part=seq[int(start/ws):int(end/ws)]
    #print("DSEQ ES", part)    
    return jsonify(response=list(part))


#%%
# positions is a dictionary of tracks -> [positions]
# types filters out annotations not corresponding to these types ("gene", "exon" etc.)
# onlyIDs if true it only returns gene ids
# align as positions are only integers, the method takes annotations in a range
#       defined as [position-window/2,position+window/2] if align="center"
#       or [position,position+window] if align="left"
# intersect if "soft" it requires that only a portion of the interval is in the annotation
#           "hard" requires the full interval to be inside the annotation 

@app.route("/annotations")
def annotations(positions={}, window=1000, types=["any"], onlyIDs="False", align="left", intersect="soft", dataName="None"):
    global data
    window=eval(request.args.get("window"))
    types=eval(request.args.get("types"))
    positions=eval(request.args.get("positions"))
    onlyIDs=str(request.args.get("onlyIDs"))
    align=str(request.args.get("align"))
    intersect=str(request.args.get("intersect"))
    dataName=str(request.args.get("dataName"))
    
    res={}
    print("---------------------------------------")
     
    gis=[]
    for track in positions.keys():
        print("Retrieving annotations on track", track)
        if(dataName=="None"):
            dbp=data["batch"]["processed"]
        else:
            dbp=data[dataName]["batch"]["processed"]
        pos=positions[track]
        if(type(pos)==str):
            pos=eval(pos)
            
        if(track in dbp["gff"].keys()):
            res[track]=annotationsLocal(positions=pos, gff=dbp["gff"][track], window=window, types=types, onlyIDs=onlyIDs, align=align, intersect=intersect)
        else:
            res[track]=[]
            
                        
    if("search" in data.keys()):
        for ch in res.keys():
            for x in res[ch].keys():
                for y in res[ch][x]:
                    if(y["t"]=="gene"):
                        gis.append(y["id"])
        data["annotations"]=res    
        data["gis"]=gis                       
    return jsonify(response=res)
#%%
def annotationsLocal(positions, gff, window=1000, types=["any"], track="None", onlyIDs="False", align="left", intersect="soft"):
    import time
    t0=time.clock()
    print("Annotating ",len(positions)," positions")
    res=ann.annotate(positions, gff, types, window, align, intersect)
    print("Annotations found for ",len(res)," positions")
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
def nucleotides(start=0, end=10, track="None", dataName="None"):
    global data

    dataName=str(request.args.get("dataName"))
    if(dataName=="None"):
        dbp=data["batch"]["processed"]
    else:
        dbp=data[dataName]["batch"]["processed"]

    start=int(request.args.get("start"))
    end=int(request.args.get("end"))
    track=str(request.args.get("track"))
    return jsonify(response=dbp["fasta"][track][start:end])
 
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
def nucProfile(positions=[], size=10, track="None", k=6, dataName="None"):
    global data
    import time
    t00=time.clock()
    pos=eval(request.args.get("positions"))
    size=int(request.args.get("size"))
    track=str(request.args.get("track"))
    k=int(request.args.get("k"))
    dataName=str(request.args.get("dataName"))

    if(dataName=="None"):
        dbp=data["batch"]["processed"]
    else:
        dbp=data[dataName]["batch"]["processed"]
    
    #basic operations (prof and c might not be very useful)
    seqs={}
    for p in pos:
        seqs[p]=dbp["fasta"][track][p:p+size].upper()
    c=ms.consensus(seqs.values())
    prof=ms.profile(seqs.values())

    #gibbs sampling for motif search
    if(len(seqs)>1):
        t0=time.clock()
        for x in range(k,k+1,1): #for different k
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
    
    dg=helpers.getDataAnnot(data)
            
    a=eval(request.args.get("genes"))
    print('annotation size: ',len(a))
    agon=ann.annotateGOnames(a, dg["goa"], dg["go"])
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

    if(len(gis)==0):
        gis=set(data["gis"])
    
    dg=helpers.getDataAnnot(data)
    
    ego={}
    ego=ann.enrichmentFisher(gis, dg["goa"], alpha, correction)
    for k in ego.keys():
        if(dg["go"].has_key(k) and dg["go"][k]!="undefined"):
            ego[k]["go_name"]=dg["go"][k]
    data["ego"]=ego
    return jsonify(response=ego)

#%% automatic assignment of users (no login/pass protocol)
@app.route("/assignUser")
def assignUser():
    import random
    import sys
    
    global user
    
    user=str(random.randrange(0, sys.maxint))
    session[user]={}
    return jsonify(response=user)
    
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
    #DISCARDED: check password and so on.
    if user in session.keys():
        data=session[user]

#@app.after_request
#def serialize_passport(response):
#    if hasattr(g, "passport"):
#        session["passport_id"] = g.passport.id
#    return response
    
@app.route("/testArray")
def testArray():
    import numpy as np
    return jsonify(a=np.abs([-0.80176, -0.098632812, 0.09375, 0.8017578125]))

@app.route("/test")
def test():
    return "Seqview server correclty configured"

#%%
@app.route("/availableOrganisms")
def availableOrganisms():
    return jsonify(response=ann.getAnnotationFolders()) 
    
#-------------------- LAUNCH -----------------
if __name__ == '__main__':
    #app.run(debug=True)
    app.run(debug=True, host='0.0.0.0', port=2750)
    
    
    