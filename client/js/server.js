/*
 ┌────────────────────────────────────────────────────────────┐
 │ server.js                                                  │
 ├────────────────────────────────────────────────────────────┤
 │ Description:                                               │
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
         * @param user
         * @param password
         * @param serverPath
         */
        Server.connect = function (DEBUG, callback, user, password, serverPath)
        {
            serverPath     || ( serverPath = "http://127.0.0.1:2750/" );


            _serverPath     = removeLastSlash(serverPath);
            _user           = user;
            _password       = password;
            _DEBUG          = DEBUG;


            // Here you should check the credentials...
            if(user == "jpiriz")
                callback(true);
            else
                callback(false);
        };


        /**
         * Uploads the file to the server. Posterior calls to REST will be by the path of the file
         * @param callback
         * @param file
         * @param forceReload
         */
        // NOTE: Remember REST Flask calls require enable CORS in the browser!!!
        Server.sendFile = function (callback, file, forceReload)
        {
            var startTime = new Date();

            // Show the image of "loading..."
            showImageLoading("imgLoadingFile", true);


            var _forceReload = "";
            if(forceReload) _forceReload = "True";
            else            _forceReload = "False";

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/testUpload?user="+_user+"&password="+_password+"&filename="+file.name+"&forceReload="+_forceReload,
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
                    if(_DEBUG) console.log("sendFile(): uploaded? "+result.response);

                    if(!(result.response == "outdated version" || result.response == "not found"))
                    {
                        // Hide the image of "loading..."
                        showImageLoading("imgLoadingFile", false);


                        if (_DEBUG) console.log("Time spent sending: " + (new Date() - startTime) + "ms");

                        callback();
                    }
                    // Only in this case, we upload it
                    else
                    {
                        // The FormData object lets you compile a set of key/value pairs to send using XMLHttpRequest (AJAX)
                        var fd = new FormData();
                        fd.append('file',file);

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
                                // Hide the image of "loading..."
                                showImageLoading("imgLoadingFile", false);

                                if (_DEBUG) console.log("sendFile(): uploaded");
                                if (_DEBUG) console.log("Time spent sending: " + (new Date() - startTime) + "ms");

                                callback();
                            })
                            .fail(function()
                            {
                                if(_DEBUG) console.log("sendFile(): upload failed...");
                            });
                    }
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("sendFile(): testUpload failed...");
                    javascript_abort();
                });
        };


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
        Server.preprocess = function (callback, filename, track, ws, nb, maxSize)
        {
            var stdev = 3;
            var recharge = "False";


            var startTime = new Date();

            // Show the image of "loading..."
            showImageLoading("imgLoadingFile", true);

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/preprocess?user="+_user+"&password="+_password+"&filename="+filename+"&track="+track+"&windowSize="+ws+"&numBins="+nb+"&maxSize="+maxSize+"&stdev="+stdev+"&recharge="+recharge,
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

                    // Hide the image of "loading..."
                    showImageLoading("imgLoadingFile", false);

                    if(_DEBUG) console.log("preprocress(): discretization done...");
                    if (_DEBUG) console.log("Time spent preprocessing: " + (new Date() - startTime) + "ms");

                    callback(response, track);
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("preprocress(): discretization failed...");
                });
        };


        /**
         * Calls the REST service available to a data summary and characterization of a particular chromosome
         * @param callback
         * @param track
         *
         */
        Server.getTrack = function (callback, track)
        {
            var startTime = new Date();

            // Show the image of "loading..."
            showImageLoading("imgLoadingFile", true);

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/getTrack?user="+_user+"&password="+_password+"&track="+track,
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

                    // Hide the image of "loading..."
                    showImageLoading("imgLoadingFile", false);

                    if(_DEBUG) console.log("getTrack(): get track done...");
                    if(_DEBUG) console.log("Time spent get track: " + (new Date() - startTime) + "ms");

                    callback(response, track);
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("getTrack(): get track failed...");
                });
        };


        /**
         * Search this pattern in all the genome.
         * @param callback
         * @param pattern
         * @param d
         */
        Server.search = function (callback, pattern, d, geo, intersect, softMutations)
        {
            var startTime = new Date();

            // Show the image of "loading..."
            showImageLoading("imgLoadingSearch", true);


            if(_DEBUG) console.log("pattern is "+pattern);
            pattern = encodeURIComponent(pattern);  // it encodes the following characters: , / ? : @ & = + $ #

            var requestAJAX = $.ajax(
                {
                    //url: _serverPath+"/search?user="+_user+"&password="+_password+"&pattern="+pattern+"&d="+d+"&geo="+geo+"&intersect="+intersect+"&softMutations="+softMutations,
                    url: _serverPath+"/search?user="+_user+"&password="+_password+"&pattern="+pattern+"&d="+d+"&geo="+geo+"&intersect=false&softMutations="+softMutations,
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
                        var response = [];
                        response.points = result.points; // points found in all chromosomes
                        response.sizePattern = result.sizePattern;

                        // Hide the image of "loading..."
                        showImageLoading("imgLoadingSearch", false);

                        if(_DEBUG) console.log("search(): search done...");
                        if(_DEBUG) console.log("Time spent searching: "+ (new Date()-startTime)+"ms");

                        callback(response);
                    }
                    else
                    {
                        if(_DEBUG) console.log("search(): "+result.msg);
                        $("#searchImg")[0].style.visibility="hidden";
                        alert("Error: "+result.msg);
                    }
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("search(): search failed...");
                    $("#searchImg")[0].style.visibility="hidden";
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
                                               gis, annotations, numChromosome, startTime, numMatches)
        {
            // Initialize variables first
            if(typeof(gis) === 'undefined')             gis = "";
            if(typeof(annotations) === 'undefined')     annotations = {};
            if(typeof(numChromosome) === 'undefined')   numChromosome = 1;
            if(typeof(startTime) === 'undefined')       startTime = new Date();
            if(typeof(numMatches) === 'undefined')      numMatches = 0;

            // Calculate points of this track and number of matches
            var track  = chromosomes[numChromosome-1];
            var points = JSON.parse(allPoints[track]);
            //for( var j=0; j<points.length; j++)
            //    points[j]*=ws;
            numMatches += points.length;
            if(_DEBUG) console.log(points.length+" matches in track: "+track);

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/annotations?user="+_user+"&password="+_password+
                    "&positions=["+points+"]&types="+types+"&window="+window+"&align="+align+"&track="+track+"&onlyIDs="+onlyIDs+"&intersect="+intersect,
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
                    var startClient=new Date()
                    console.log("Received at: "+startClient);
                    numChromosome++;
                    if(gis != "" && result.response != "") gis += ",";

                    if(onlyIDs=="True")
                        gis += result.response+",";
                    else
                    {
                        var response = result.response;
                        for (var key in response)
                        {
                            // for...in loops over all enumerable properties
                            // (which is not the same as "all array elements"!)
                            if(response.hasOwnProperty(key))
                            {
                                var rrk = result.response[key];
                                for (var i in rrk)
                                {
                                    if(rrk.hasOwnProperty(i))
                                    {
                                        var rrki = result.response[key][i];
                                        gis += rrki["id"]+",";
                                        annotations[rrki["id"]] = {};
                                        annotations[rrki["id"]]["pos"] = key;
                                        annotations[rrki["id"]]["name"] = rrki["n"];
                                        annotations[rrki["id"]]["sense"] = rrki["ss"];
                                        annotations[rrki["id"]]["start"] = rrki["s"];
                                        annotations[rrki["id"]]["end"] = rrki["e"];
                                        annotations[rrki["id"]]["chromosome"] = track;
                                    }
                                }
                            }
                        }
                    }
                    console.log("Time spent in data processing "+(new Date()-startClient)+"ms")

                    // While the number of chromosomes is less than or equal, we continue to get more annotations
                    if(numChromosome <= chromosomes.length)
                    {

                        var startWS=new Date()
                        Server.allAnnotationsGenes(getEnrichment, allPoints, "[\"gene\"]", window, "left", "False", chromosomes, GVB_GLOBAL.ws, intersect,
                                                    gis, annotations, numChromosome, (new Date()-startTime), numMatches);
                        console.log("Time spent getting annotations: "+ (new Date()-startWS)+"ms");
                    }
                    else
                    {
                        if(_DEBUG) console.log("Total number of matches: "+numMatches);
                        if(_DEBUG) console.log("allAnnotationsGenes(): get all annotations done...");
                        if(_DEBUG) console.log("Time spent getting all annotations: "+ (new Date()-startTime)+"ms");

                        callback(gis, annotations);
                    }
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("allAnnotationsGenes(): get all annotations failed...");
                });
        };


        /**
         * Get enrichment in all the genome (with 'gis').
         * @param callback
         * @param gis
         * @param correction  ('none', 'bonferroni', 'fdr', 'fwer')
         * @param alpha
         */
        Server.enrichmentGO = function (callback, gis, correction, alpha)
        {
            var startTime = new Date();

            //Not sure if this duplicate removal for sending is really healping or not (seems not)
            var arrayGis0 = gis;
            arrayGis0 = arrayGis0.split(',');
            var arrayGis=[]
            for(var i in arrayGis0)
                {
                if(arrayGis.indexOf(arrayGis0[i])<0)
                    arrayGis.push(arrayGis0[i]);
                }
            console.log("Number of gene annotations: "+arrayGis0.length+", of which unique: "+arrayGis.length);
            gis=arrayGis.join(",");
            //

            // Parse 'gis' to convert it to array: a,b,c,d => ["a","b","c","d"]
            gis = "[\""+gis+"\"]";
            gis = gis.replace(/,/g,"\",\"");  //replaces ',' by '","'


            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/enrichmentGO?user="+_user+"&password="+_password+"&annotations="+gis+"&correction="+correction+"&alpha="+alpha,
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
                .fail(function()
                {
                    if(_DEBUG) console.log("enrichmentGO(): enrichmentGO failed...");
                });
        };


        /**
         * Get genes of interest of a particular chromosome
         * @param callback
         * @param point
         * @param types
         * @param window
         * @param align
         * @param track
         * @param onlyIDs
         */
        Server.annotationsGenes = function (callback, point, types, window, align, track, onlyIDs)
        {
            // var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/annotations?user="+_user+"&password="+_password+
                                    "&positions=["+point+"]&types="+types+"&window="+window+"&align=\""+align+"\"&track="+track+"&onlyIDs="+onlyIDs+"&intersect="+GVB_GLOBAL.intersect,
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
                    var annotations = result.response;

                    // console.log("annotationsGenes(): get of gene done...");
                    // console.log("Time spent annotationsGenes: "+ (new Date()-startTime)+"ms");

                    if(annotations.hasOwnProperty(point)) // if object "annotations" has no the point, don't draw the line annotations
                    {
                        callback(annotations[point]);
                    }
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("annotationsGenes(): annotationsGenes failed...");
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
                    if(false) console.log("Time spent annotationsGenes: "+ (new Date()-startTime)+"ms");
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
                .fail(function()
                {
                    if(_DEBUG) console.log("annotationsGenes(): annotationsGenes failed...");
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
        Server.getPartSeq = function (callback, track, startSeq, endSeq, numNucleotides, point, sizePattern)
        {
            // var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/getPartSeq?user="+_user+"&password="+_password+"&track="+track+"&start="+startSeq+"&end="+endSeq,
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
                    response.partSeq = result.partSeq;

                    callback(response.partSeq, numNucleotides, point, sizePattern);
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("getPartSeq(): getPartSeq failed...");
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
        Server.nucleotides = function (callback, track, startSeq, endSeq, start, point)
        {
            // var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/nucleotides?user="+_user+"&password="+_password+"&track="+track+"&start="+startSeq+"&end="+endSeq,
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

                    callback(start, point, response)
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("nucleotides(): get nucleotides failed...");
                });
        };


        /**
         * Calls the REST service available to retrieve a sequence of nucleotides
         * @param callback
         * @param positions: array with numerical starting positions
         * @param size: length from the starting points
         * @param track: chromosome or track where the sequences must be taken from
         */
        Server.nucProfile = function (callback, track, positions, size)
        {
            // var startTime = new Date();

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/nucProfile?user="+_user+"&password="+_password+"&track="+track+"&positions="+positions+"&size="+size,
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
                .fail(function()
                {
                    if(_DEBUG) console.log("nucProfile(): discretization failed...");
                });
        };


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
            msg   || ( msg = 'This is not an error. This is just to abort javascript' );
            alert(msg);

            // Hide all images of "loading..."
            showImageLoading("imgLoadingFile", false);
            showImageLoading("imgLoadingSearch", false);

            throw "ERROR GBV: "+msg;
        }

        return Server;
    }

})(window);



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
