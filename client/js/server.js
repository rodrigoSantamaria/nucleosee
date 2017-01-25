/*
 ┌────────────────────────────────────────────────────────────┐
 │ server.js                                                  │
 ├────────────────────────────────────────────────────────────┤
 │ Description: Model logic for the browser (backend communication)
 |  GPL v3 by Rodrigo Santamaría (University of Salamanca)    |
 └────────────────────────────────────────────────────────────┘
 */


// To try calls:
/*
curl -i -H "Accept: application/json" -H "Content-Typ: application/json" -X GET \
"http://127.0.0.1:2750/testUpload?user=jpiriz&password=ninguna&filename=23479_h90_wlt_mean.wig&forceReload=false"
*/


(function(window)
{
    // I recommend use this:
    // http://stackoverflow.com/questions/1335851/what-does-use-strict-do-in-javascript-and-what-is-the-reasoning-behind-it
    'use strict';

    // Define globally if it doesn't already exist
    if(typeof(Server) === 'undefined')
    {
        window.Server = defineServerLibrary();
    }
    else
    {
        console.log("Library \"Server\" already defined.");
    }


    function defineServerLibrary()
    {
        var Server = {};

        var _DEBUG      = false;
        var _serverPath = "";
        var _user       = "";
        var _password   = "";


        /**
         * Connect with the server.
         * @param DEBUG
         * @param callback
         * @param user      User and password are deprecated
         * @param password
         * @param serverPath    Change to your server on index.html call
         */
        Server.connect = function (DEBUG, callback, user, password, serverPath)
        {
            serverPath     || ( serverPath = "http://127.0.0.1:2750/" );


            _serverPath     = removeLastSlash(serverPath);
            _user           = user;
            _password       = password;
            _DEBUG          = DEBUG;

            // We do just a simple user random assigning to control concurrent accesses
            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/assignUser",
                    type: "GET",
                    datatype: "json"
                });


            $.when(requestAJAX)
                .done(function(result)
                {
                    _user=result.response;
                    callback(true);
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("connect(): Connect failed...");
                    javascript_abort("Connection failed. Server might be down: "+errorThrown+" "+jqXHR.responseText);
                    callback(false);
                });
        };


        /**
         * Uploads the file to the server. Posterior calls to REST will be by the path of the file
         * @param file
         * @param forceReload
         */
        Server.testFile = function (file)
        {
            var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    //url: _serverPath+"/testUpload?user="+_user+"&password="+_password+"&filename="+file.name+"&forceReload="+_forceReload,
                    url: _serverPath+"/testUpload?user="+_user+"&password="+_password+"&filename="+encodeURIComponent(file.name)+"&forceReload=false",
                    type: "GET",
                    datatype: "json"
                });

            //noinspection JSUnresolvedFunction
            /**
             * @typedef {Object} result
             * @property response
             */
            $.when(requestAJAX)
                .done(function(result)
                {
                    if(_DEBUG) console.log("testFile(): uploaded? "+result.response);
                    GVB_GLOBAL.sendFile=false;
                    if(!(result.response == "outdated version" || result.response == "not found"))
                        {
                        // Hide the image of "loading..."
                        showImageLoading("imgLoadingFile", false);

                        if (_DEBUG) console.log("Time spent sending: " + (new Date() - startTime) + "ms");

                        //callback();
                        //Change load options if there's already a preprocessing done
                        $("#paramSD")[0].value=result.clip;
                        $("#paramWS")[0].value=result.ws;
                        $("#paramNB")[0].value=result.nb;
                        $("#speciesList")[0].value=result.org;
                        $("#interpolationList")[0].value=result.interpol.toLowerCase();

                        $("#btLoadOptions")[0].disabled=false;
                        }
                    if(result.response=="not found")
                        GVB_GLOBAL.sendFile=true;
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("testFile(): testUpload failed...");
                    console.log(jqXHR)
                    javascript_abort("Error checking the file. "+errorThrown+" "+jqXHR.responseText);
                });
        };



        Server.sendFile=function(callback, files, i)
            {
            // The FormData object lets you compile a set of key/value pairs to send using XMLHttpRequest (AJAX)
            var fd = new FormData();
            fd.append('file',files[i]);

            var requestAJAX2 = $.ajax(
                {
                    url: _serverPath+"/upload?user="+_user+"&password="+_password,
                    type: "POST",
                    data: fd,
                    processData: false, // tell jQuery not to process the data
                    contentType: false  // tell jQuery not to set contentType
                });

            //noinspection JSUnresolvedFunction
            $.when(requestAJAX2)
                .done(function()
                {

                    if (_DEBUG) console.log("testFile(): uploaded");

                    if(i<files.length-1)
                        Server.sendFile(callback, files, i+1)
                    else {
                        // Hide the image of "loading..."
                        showImageLoading("imgLoadingFile", false);
                        callback();
                        }
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("testFile(): upload failed...");
                    console.log(jqXHR)
                    javascript_abort("Error uploading the file. Available formats are .wig and .bw "+errorThrown+" "+jqXHR.responseText);
                });
            }


        /**
         * Calls the REST service available to a summary and characterization of data
         * @param callback
         * @param filename  usually a .wig file
         * @param track
         * @param ws    window size for summarization (usually 100 or more)
         * @param nb    number of bins for summarization (typically 5 or 7)
         * @param maxSize maximum size of the seq array that we want
         * @returns {Array} with the following fields
         *          seq - original (normalized) sequence, sampled to fit into maxSize array
         *          max - maximum value in seq
         *          min - minimum value in seq
         *          mean - mean value in seq
         *          stdev -  standard deviation in seq
         *          dseq - each ws values in seq are summarized as a letter depending of its average
         *                  e.g. if we have a seq with mean=10 and sd=5 and we have a ws=3, and nb=5 a window of values
         *                      10,9,10 will be coded as 'c' and another one with values 16,20,12 will be coded as 'e' and so
         *                      on.
         */
        Server.preprocess = function (callback, filenames, track, ws, nb, maxSize, organism, interpolation, stdev, description)
        {
            var startTime = new Date();

            // Show the image of "loading..."
            showImageLoading("imgLoadingFile", true);

            for(var i in filenames)
                filenames[i] = "\""+encodeURIComponent(filenames[i])+"\"";  // it encodes the following characters: , / ? : @ & = + $ #

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/preprocess?user="+_user+"&password="+_password+"&filenames=["+filenames+"]&track="+track+"&windowSize="+ws+"&numBins="+nb+"&maxSize="+maxSize+"&stdev="+stdev+"&dataName="+description+"&organism="+organism+"&interpolation="+interpolation,
                    type: "GET",
                    datatype: "json"
                });

            //noinspection JSUnresolvedFunction
            /**
             * @typedef {Object} result
             * @property seq
             * @property maximum
             * @property minimum
             * @property mean
             * @property stdev
             * @property dseq
             * @property fullLength
             * @property chromosomes
             */
            $.when(requestAJAX)
                .done(function(result)
                {
                    var response = [];
                    response.seq=result.seq; // this is only a sample, as it is too large to show as a whole and to send via REST
                    response.max=result.maximum;
                    response.min=result.minimum;
                    response.mean=result.mean;
                    response.stdev=result.stdev;
                    response.dseq=result.dseq;
                    response.fullLength=result.fullLength;
                    response.chromosomes=result.chromosomes;
                    response.bins=result.bins;
                    response.filenames=filenames;

                    // Hide the image of "loading..."
                    showImageLoading("imgLoadingFile", false);

                    //This is to force updating of data list
                    document.getElementById('selectionList').options.length=0;


                    if(_DEBUG) console.log("preprocress(): discretization done...");
                    if (_DEBUG) console.log("Time spent preprocessing: " + (new Date() - startTime) + "ms");

                    callback(response, track,0,1);
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("preprocress(): discretization failed...");
                    javascript_abort("preprocessing() failed: possible causes: \n· Make sure the wig file follows the standard format\n· Make sure you selected the right species");
                });
        };


        /**
         * Lists available data in the server
         */
        Server.listData = function (callback)
        {
            var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/listData?user="+_user+"&password="+_password,
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                    var response = result.response;
                    if(result.response=="error")
                        {
                        javascript_abort(result.msg);
                        return;
                        }


                    if(_DEBUG) console.log("preprocress(): discretization done...");
                    if (_DEBUG) console.log("Time spent preprocessing: " + (new Date() - startTime) + "ms");

                    callback(response);
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("listData(): listing failed...");
                    javascript_abort("listData() failed: possibly due to corruption in the server, ask your administrator");
                });
        };


        /**
         * Removes user from server
         */
        Server.removeUser = function ()
        {

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/removeUser?user="+_user+"&password="+_password,
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                   console.log("User succesfully removed from the server")
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    console.log("removeUser() failed...");
                });
        };


        Server.getDSeq = function (start,end,track,dataName)
        {
            var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/getDSeq?user="+_user+"&password="+_password+"&start="+start+"&end="+end
                        +"&dataName="+dataName+"&track="+track,
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                    var response = result.response;
                    console.log("DSEQ("+dataName+")"+response)
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("getDSeq(): listing failed...");
                    javascript_abort("getDSeq() failed: possibly due to corruption in the server, ask your administrator");
                });
        };

        /**
         * Selects available preprocessed data in the server
         */
        Server.selectData = function (callback, dataName, index, total, clear)
        {
            var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/loadData?user="+_user+"&password="+_password+"&dataName="+dataName+"&clear="+clear,
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                    var response = [];
                    response.seq=result.seq; // this is only a sample, as it is too large to show as a whole and to send via REST
                    response.max=result.maximum;
                    response.min=result.minimum;
                    response.mean=result.mean;
                    response.stdev=result.stdev;
                    response.dseq=result.dseq;
                    response.fullLength=result.fullLength;
                    response.chromosomes=result.chromosomes;
                    response.bins=result.bins;
                    response.dataName=dataName;
                    GVB_GLOBAL.ws=result.windowSize;

                    if(_DEBUG) console.log("preprocress(): discretization done...");
                    if (_DEBUG) console.log("Time spent preprocessing: " + (new Date() - startTime) + "ms");

                    callback(response, "None", index, total);
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("listData(): listing failed...");
                    javascript_abort("listData() failed: possibly due to an error in the server, ask your administrator");
                });
        };

        /**
         * Calls the REST service available to a data summary and characterization of a particular chromosome
         * @param callback
         * @param track
         *
         */
        Server.getTrack = function (callback, track, dataName, index, total)
        {
            var startTime = new Date();

            // Show the image of "loading..."
            showImageLoading("imgLoadingFile", true);

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/getTrack?user="+_user+"&password="+_password+"&track="+track+"&dataName="+dataName,
                    type: "GET",
                    datatype: "json"
                });

            //noinspection JSUnresolvedFunction
            /**
             * @typedef {Object} result
             * @property seq
             * @property maximum
             * @property mean
             * @property stdev
             * @property dseq
             * @property fullLength
             * @property chromosomes
             * @property search
             */
            $.when(requestAJAX)
                .done(function(result)
                {
                    var response = [];
                    response.seq=result.seq; // this is only a sample, as it is too large to show as a whole and to send via REST
                    response.max=result.maximum;
                    response.min=result.minimum;
                    response.mean=result.mean;
                    response.stdev=result.stdev;
                    response.dseq=result.dseq;
                    response.fullLength=result.fullLength;
                    response.bins=result.bins;
                    response.chromosomes=result.chromosomes;
                    response.search=result.search;
                    response.ego=result.ego;
                    response.dataName=dataName;

                    // Hide the image of "loading..."
                    showImageLoading("imgLoadingFile", false);

                    if(_DEBUG) console.log("getTrack(): get track done...");
                    if(_DEBUG) console.log("Time spent get track: " + (new Date() - startTime) + "ms");

                    callback(response, track, index, total);

                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("getTrack(): get track failed...");
                    javascript_abort("getTrack() failed");

                });
        };


        /**
         * Search this pattern in all the genome.
         * @param callback
         * @param pattern
         * @param d
         */
        Server.search = function (callback, pattern, d, geo, intersect, softMutations, dataName1,dataName2,join,pattern2)
        {
            var startTime = new Date();

            // Show the image of "loading..."
            showImageLoading("imgLoadingSearch", true);


            if(_DEBUG) console.log("pattern is "+pattern);
            pattern = encodeURIComponent(pattern);  // it encodes the following characters: , / ? : @ & = + $ #

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/search?user="+_user+"&password="+_password+"&pattern="+pattern+"&d="+d+
                        "&geo="+geo+"&intersect="+intersect+"&softMutations="+softMutations+
                        "&dataName1="+dataName1+"&dataName2="+dataName2+"&join="+join+"&pattern2="+pattern2,
                    type: "GET",
                    datatype: "json"
                });

            //noinspection JSUnresolvedFunction
            /**
             * @typedef {Object} result
             * @property points
             * @property sizePattern
             * @property msg
             */
            $.when(requestAJAX)
                .done(function(result)
                {
                    if(result.response != "error")
                    {
                        /*var response = [];
                        response.points = result.points; // points found in all chromosomes
                        response.sizePattern = result.sizePattern;*/
                        var response=result;

                        // Hide the image of "loading..."
                        showImageLoading("imgLoadingSearch", false);

                        if(_DEBUG) console.log("search(): search done...");
                        if(_DEBUG) console.log("Time spent searching: "+ (new Date()-startTime)+"ms");

                        callback(response);
                    }
                    else
                    {
                        if(_DEBUG) console.log("search(): "+result.msg);
                        javascript_abort("Search failed: "+result.msg);
                    }
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("search(): search failed...");
                    javascript_abort("search() failed");
                });
        };


        /**
         * Get genes of interest in all the genome.
         * @param callback
         * @param allPoints
         * @param types
         * @param window
         * @param align
         * @param onlyIDs
         * @param chromosomes
         * @param ws
         * @param gis               it's not necessary!! sure?
         * @param annotations       it's not necessary!! sure?
         * @param numChromosome     it's not necessary!!
         * @param startTime         it's not necessary!!
         * @param numMatches        it's not necessary!!
         */
        Server.allAnnotationsGenes = function (callback, allPoints, types, window, align, onlyIDs, chromosomes, ws, intersect,
                                               dataName, track)
        {
            // Initialize variables first
            var gis = "";
            var annotations = {};
            var startTime = new Date();

            // Calculate points of this track and number of matches
            //var points = JSON.parse(allPoints);
            var points = JSON.stringify(allPoints);
            var ws=window;
            if(typeof(window)=="object")
                ws=JSON.stringify(window);

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/annotations?user="+_user+"&password="+_password+
                    "&positions="+points+"&types="+types+"&window="+ws+"&align="+align+
                    "&onlyIDs="+onlyIDs+"&intersect="+intersect+"&dataName="+dataName,
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                    var startClient=new Date()
                    console.log("Received at: "+startClient);
                    if(gis != "" && result.response != "") gis += ",";

                    if(onlyIDs=="True")
                        gis += result.response+",";
                    else
                    {
                        var response = result.response;
                        for (var key in response)//for each chromosome or track
                        {
                            // for...in loops over all enumerable properties
                            // (which is not the same as "all array elements"!)
                            if(response.hasOwnProperty(key))
                            {
                                var rrk = result.response[key];
                                for (var i in rrk)  //for each position with annotations
                                {
                                    if(rrk.hasOwnProperty(i))
                                    {
                                        var rrki = result.response[key][i];
                                        for(var k in rrki) {    //for each annotation at a given position
                                            var rrkik=rrki[k]
                                            if(rrkik["t"]=="gene") {
                                                var id=rrkik["id"]
                                                gis += id + ",";
                                                annotations[id] = {};
                                                annotations[id]["pos"] = i;
                                                annotations[id]["name"] = rrkik["n"];
                                                annotations[id]["sense"] = rrkik["ss"];
                                                annotations[id]["start"] = rrkik["s"];
                                                annotations[id]["end"] = rrkik["e"];
                                                annotations[id]["chromosome"] = key;
                                                }
                                            }
                                    }
                                }
                            }
                        }
                    }
                    // While the number of chromosomes is less than or equal, we continue to get more annotations
                    if(_DEBUG) console.log("allAnnotationsGenes(): get all annotations done...");
                    if(_DEBUG) console.log("Time spent getting all annotations: "+ (new Date()-startTime)+"ms");

                    callback(gis, annotations);

                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("allAnnotationsGenes(): get all annotations failed...");
                    javascript_abort("allAnnotationsGenes failed");

                });
        };



        /**
         * Get enrichment in all the genome (with 'gis').
         * @param callback
         * @param gis: leave empty ("") as it is now stored in the backend
         * @param correction  ('none', 'bonferroni', 'fdr', 'fwer')
         * @param alpha
         */
        Server.enrichmentGO = function (callback, gis, correction, alpha, discard)
        {
            var startTime = new Date();

            if(gis.length>0)
                {
                //Not sure if this duplicate removal for sending is really healping or not (seems not)
                var arrayGis0 = gis;
                arrayGis0 = arrayGis0.split(',');
                var arrayGis = []
                for (var i in arrayGis0) {
                    if (arrayGis.indexOf(arrayGis0[i]) < 0)
                        arrayGis.push(arrayGis0[i]);
                }
                console.log("Number of gene annotations: " + arrayGis0.length + ", of which unique: " + arrayGis.length);
                gis = arrayGis.join(",");
                //

                // Parse 'gis' to convert it to array: a,b,c,d => ["a","b","c","d"]
                gis = "[\"" + gis + "\"]";
                gis = gis.replace(/,/g, "\",\"");  //replaces ',' by '","'
                }
            else
                gis="[]"


            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/enrichmentGO?user="+_user+"&password="+_password+"&annotations="+gis+"&correction="+correction+"&alpha="+alpha+"&discard="+discard,
                    type: "GET",
                    datatype: "json"
                });

            //noinspection JSUnresolvedFunction
            /**
             * @typedef {Object} result
             * @property response
             */
            $.when(requestAJAX)
                .done(function(result)
                {
                    var response = result.response;

                    if(_DEBUG) console.log("enrichmentGO(): enrichmentGO done...");
                    if(_DEBUG) console.log("Time spent Enrichment analysis: "+ (new Date()-startTime)+"ms");

                    callback(response);

                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("enrichmentGO(): enrichmentGO failed...");
                    javascript_abort("enrichmentGO() failed, check that python fisher library is installed in your server");

                });
        };


        /**
         * Get the annotations for a single point
         * @param callback
         * @param point Must be referenced as a dictionary chromosome->[point]
         * @param types
         * @param window
         * @param align
         * @param track
         * @param onlyIDs
         */
        Server.annotationsPoint = function (callback, point, types, window, align, onlyIDs, dataName)
        {
            var positions=JSON.stringify(point)
            var p;
            for(var ch in point)
                for(var i in point[ch])
                    {
                    p = point[ch][i]
                    break;
                    }

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/annotations?user="+_user+"&password="+_password+
                                    "&positions="+positions+"&types="+types+"&window="+window+"&align=\""+align+
                                    "&onlyIDs="+onlyIDs+"&intersect="+GVB_GLOBAL.intersect+
                                    "&dataName="+dataName,
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                    var annotations = result.response;

                    for(var ch in annotations)
                        {
                        if (annotations[ch].hasOwnProperty(p)) // if object "annotations" has no the point, don't draw the line annotations
                            callback(annotations[ch][p]);
                        }
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("annotationsPoint(): annotationsPoint failed...");
                    javascript_abort("annotationsPoint() failed");

                });
        };


        /**
         * Get GO terms for a given list of gene ids
         * @param callback
         * @param genes list of genes
         */
        Server.Gene2GO = function (callback, genes)
        {
            var startTime = new Date();
            var genes1 = "[\""+genes+"\"]";
            genes1 = genes1.replace(/,/g,"\",\"");  //replaces ',' by '","'

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/annotationsGOA?user="+_user+"&password="+_password+
                    "&genes="+genes1+"&types=[\"any\"]",
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                    var annotations = result.response;

                    if(false) console.log("Gene2GO(): get of gene done...");
                    if(false) console.log("Time spent annotationsPoint: "+ (new Date()-startTime)+"ms");
                    var goterms={};
                    for(var i in genes) {
                        var g = genes[i];
                        goterms[g] = [];
                    }
                    for (var key in annotations)
                        {
                        for(var i in annotations[key]["genes"])
                            goterms[annotations[key]["genes"][i]].push(annotations[key]["name"]);
                        }
                    callback(genes, goterms);
                    //return annotations;
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("Gene2GO(): failed...");
                    javascript_abort("Gene2GO() failed");

                });
        };


        /**
         * Calls the REST service available to a summary and characterization of data (more reduced)
         * @param callback
         * @param track
         * @param startSeq
         * @param endSeq
         * @param numNucleotides
         * @param point
         * @param sizePattern
         */
        Server.getPartSeq = function (callback, track, startSeq, endSeq, numNucleotides, point, sizePattern, dataName)
        {
            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/getPartSeq?user="+_user+"&password="+_password+"&track="+track+
                        "&start="+startSeq+"&end="+endSeq+"&dataName="+dataName,
                    type: "GET",
                    datatype: "json"
                });

            //noinspection JSUnresolvedFunction
            /**
             * @typedef {Object} result
             * @property partSeq
             */
            $.when(requestAJAX)
                .done(function(result)
                {
                    var response = [];
                    response=result;
                    callback(response, numNucleotides, point, sizePattern, dataName);
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("getPartSeq(): getPartSeq failed...");
                    javascript_abort("getPartSeq() failed");

                });
        };

        Server.getAllPartSeq = function (callback, start, end, point, dl1, dl2, globalSeqs, i, seqs)
        {
            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/getPartSeq?user="+_user+"&password="+_password+"&track="+globalSeqs[i].track+
                    "&start="+start+"&end="+end+"&dataName="+globalSeqs[i].dataName,
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                var response = {};
                response = result;
                seqs[globalSeqs[i].dataName]=response;


                if(i>=globalSeqs.length-1) {
                    callback(seqs, start, end);//TODO: implement callback on draw.js
                    }
                else
                    {
                    Server.getAllPartSeq(callback, start, end, point, dl1, dl2, globalSeqs, i+1, seqs)
                    }
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("getPartSeq(): getPartSeq failed...");
                    javascript_abort("getPartSeq() failed");

                });
        };

        /**
         * Gets all available organisms in the server
         * @param callback
         */
        Server.availableOrganisms = function (callback) {
        var requestAJAX = $.ajax(
            {
                url: _serverPath+"/availableOrganisms?user="+_user+"&password="+_password,
                type: "GET",
                datatype: "json"
            });

        /**
         * @typedef {Object} result
         * @property response
         */
        $.when(requestAJAX)
            .done(function(result)
            {
                callback(result.response)
            })
            .fail(function(jqXHR, textStatus, errorThrown)
            {
                if(_DEBUG) console.log("availableOrganisms(): availableOrganisms failed...");
                javascript_abort("availableOrganisms() failed");
            });

    };

            /**
         * Calls the REST service available to retrieve a sequence of nucleotides
         * @param callback
         * @param track
         * @param startSeq
         * @param endSeq
         * @param start
         * @param point
         * @returns sequence for that interval
         */
        Server.nucleotides = function (callback, track, startSeq, endSeq, start, point, dataName)
        {
            // var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/nucleotides?user="+_user+"&password="+_password+"&track="+
                    track+"&start="+startSeq+"&end="+endSeq+"&dataName="+dataName,
                    type: "GET",
                    datatype: "json"
                });

            //noinspection JSUnresolvedFunction
            /**
             * @typedef {Object} result
             * @property response
             */
            $.when(requestAJAX)
                .done(function(result)
                {
                    var response = [];
                    response.seq=result.response; // this is only a sample, as it is too large to show as a whole and to send via REST

                    // console.log("nucleotides(): get nucleotides done...");
                    // console.log("Time spent nucleotides: "+ (new Date()-startTime)+"ms");

                    callback(start, point, response, dataName)
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("nucleotides(): get nucleotides failed...");
                    javascript_abort("nucleotides() failed");

                });
        };


        /**
         * Calls the REST service available to retrieve a sequence of nucleotides
         * @param callback
         * @param positions: array with numerical starting positions
         * @param size: length from the starting points
         * @param track: chromosome or track where the sequences must be taken from
         */
        Server.nucProfile = function (callback, track, positions, size, k, dataName)
        {
            // var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/nucProfile?user="+_user+"&password="+_password+"&track="+track
                        +"&positions="+positions+"&size="+size+"&k="+k+"&dataName="+dataName,
                    type: "GET",
                    datatype: "json"
                });

            //noinspection JSUnresolvedFunction
            /**
             * @typedef {Object} result
             * @property seqs
             */
            $.when(requestAJAX)
                .done(function(result)
                {
                    var response = [];
                    response.seqs      = result.seqs;
                    response.profile   = result.profile;
                    response.consensus = result.consensus;

                    // The following three may be undefined if the number of positions is larger than 50-100
                    response.alignment  = result.alignment;
                    response.aconsensus = result.aconsensus;
                    response.aprofile   = result.aprofile;

                    // These ones are related to motif finding
                    response.motifs         = result.motifs;
                    response.locations      = result.locations;  // Motifs locations inside seqs
                    response.motifConsensus = result.motifConsensus;
                    response.motifProfile   = result.motifProfile;

                    // console.log("nucProfile(): done...");
                    // console.log("Time spent nucProfile: "+ (new Date()-startTime)+"ms");

                    callback(response);
                })
                .fail(function(jqXHR, textStatus, errorThrown)
                {
                    if(_DEBUG) console.log("nucProfile(): discretization failed...");
                    console.log(requestAJAX.url);
                    javascript_abort("nucProfile failed");
                });
        };



        function removeLastSlash(url)
        {
            var urlEnd = url.length - 1;
            if(url.lastIndexOf('/') === urlEnd) {
                url = url.substring(0, urlEnd);
            }
            return url;
        }

        function javascript_abort(msg)
        {
            msg   || ( msg = 'An error has occurred' );
            alert(msg);

            // Hide all images of "loading..."
            showImageLoading("imgLoadingFile", false);
            showImageLoading("imgLoadingSearch", false);

            throw "ERROR GBV: "+msg;
        }

        return Server;
    }

})(window);

function showImageLoading(idImageLoading, isVisible)
{
    var imageLoading = $("#"+idImageLoading);

    if(isVisible)
    {
        imageLoading.css('visibility', 'visible');
    }
    else
    {
        imageLoading.css('visibility', 'hidden');
    }
}


(function() {

    // get a reference to the d3.selection prototype,
    // and keep a reference to the old d3.selection.on
    var d3_selectionPrototype = d3.selection.prototype,
        d3_on = d3_selectionPrototype.on;

    // our shims are organized by event:
    // "desired-event": ["shimmed-event", wrapperFunction]
    var shims = {
        "mouseenter": ["mouseover", relatedTarget],
        "mouseleave": ["mouseout", relatedTarget]
    };

    // rewrite the d3.selection.on function to shim the events with wrapped
    // callbacks
    d3_selectionPrototype.on = function(evt, callback, useCapture) {
        var bits = evt.split("."),
            type = bits.shift(),
            shim = shims[type];
        if (shim) {
            evt = [shim[0], bits].join(".");
            callback = shim[1].call(null, callback);
            return d3_on.call(this, evt, callback, useCapture);
        } else {
            return d3_on.apply(this, arguments);
        }
    };

    function relatedTarget(callback) {
        return function() {
            var related = d3.event.relatedTarget;
            if (this === related || childOf(this, related)) {
                return undefined;
            }
            return callback.apply(this, arguments);
        };
    }

    function childOf(p, c) {
        if (p === c) return false;
        while (c && c !== p) c = c.parentNode;
        return c === p;
    }

})();
