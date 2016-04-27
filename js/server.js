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
"http://127.0.0.1:5000/testUpload?user=jpiriz&password=ninguna&filename=23479_h90_wlt_mean.wig&forceReload=false"
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
        console.log("Library Server already defined.");
    }


    function defineServerLibrary()
    {
        var Server = {};

        var _serverPath = "";
        var _user       = "";
        var _password   = "";
        var _DEBUG      = false;
        var _loadingImage = null;


        Server.connect = function (DEBUG, user, password, loadingImage, serverPath)
        {
            serverPath     || ( serverPath = "http://127.0.0.1:5000/" );
            loadingImage   || ( loadingImage = null );


            _serverPath     = removeLastSlash(serverPath);
            _user           = user;
            _password       = password;
            _DEBUG          = DEBUG;
            if(loadingImage != null)
            {
                _loadingImage   = $('#'+loadingImage);
            }
            else
            {
                _loadingImage = null;
            }


            // TODO: si no puede conectar, que pare todo y lance un error
            if(user == "jpiriz")
                return true;
            else
                return false;
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
            var startTime=new Date();

            // Show the image of "loading...".
            if(_loadingImage != null) { _loadingImage.css("top", 0).show(); }

            var _forceReload = "";
            if(forceReload) _forceReload = "True";

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/testUpload?user="+_user+"&password="+_password+"&filename="+file.name+"&forceReload="+_forceReload,
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                    if(_DEBUG) console.log("sendFile(): uploaded? "+result.response);

                    if(!(result.response == "outdated version" || result.response == "not found"))
                    {
                        // Hide the image of "loading...".
                        if(_loadingImage != null) { _loadingImage.css("display", "none"); }

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

                        $.when(requestAJAX2)
                            .done(function(result)
                            {
                                // Hide the image of "loading...".
                                if(_loadingImage != null) { _loadingImage.css("display", "none"); }

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
            var startTime = new Date();

            // Show the image of "loading...".
            if(_loadingImage != null) { _loadingImage.css("top", 0).show(); }

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/preprocess?user="+_user+"&password="+_password+"&filename="+filename+"&track="+track+"&windowSize="+ws+"&numBins="+nb+"&maxSize="+maxSize,
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

                    // Hide the image of "loading...".
                    if(_loadingImage != null) { _loadingImage.css("display", "none"); }

                    if(_DEBUG) console.log("preprocress(): discretization done...");
                    if (_DEBUG) console.log("Time spent preprocessing: " + (new Date() - startTime) + "ms");

                    callback(response, track);
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("preprocress(): discretization failed...");
                });
        };


        Server.getTrack = function (callback, track)
        {
            var startTime = new Date();

            // Show the image of "loading...".
            if(_loadingImage != null) { _loadingImage.css("top", 0).show(); }

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/getTrack?user="+_user+"&password="+_password+"&track="+track,
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

                    // Hide the image of "loading...".
                    if(_loadingImage != null) { _loadingImage.css("display", "none"); }

                    if(_DEBUG) console.log("getTrack(): get track done...");
                    if (_DEBUG) console.log("Time spent get track: " + (new Date() - startTime) + "ms");

                    callback(response, track);
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("getTrack(): get track failed...");
                });
        };


        Server.search = function (callback, pattern, d)
        {
            var startTime = new Date();

            // Show the image of "loading...".
            if(_loadingImage != null) { _loadingImage.css("top", 0).show(); }

            if(_DEBUG) console.log("pattern is "+pattern);
            pattern = encodeURIComponent(pattern);  // it encodes the following characters: , / ? : @ & = + $ #

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/search?user="+_user+"&password="+_password+"&pattern="+pattern+"&d="+d,
                    type: "GET",
                    datatype: "json"
                });


            $.when(requestAJAX)
                .done(function(result)
                {
                    if(result.response != "error")
                    {
                        // Convert to array...
                        var response = [];
                        response.points = result.points; // points found in all chromosomes
                        response.sizePattern = result.sizePattern;

                        // Hide the image of "loading...".
                        if(_loadingImage != null) { _loadingImage.css("display", "none"); }

                        if(_DEBUG) console.log("search(): search done...");
                        if(_DEBUG) console.log("Time spent searching: "+ (new Date()-startTime)+"ms");

                        callback(response);
                    }
                    else
                    {
                        if(_DEBUG) console.log("search(): "+result.msg);
                    }
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("search(): search failed...");
                });
        };


        Server.annotationsGenes = function (points, types, window, align, track, onlyIDs)
        {
            var response = [];
            $.ajax(
                {
                    // TODO: comprobar que align con " funciona bien y si se puede quitar
                    url: _serverPath+"/annotations?user="+_user+"&password="+_password+"&positions=["+points+"]&types="+types+"&window="+window+"&align=\""+align+"\"&track="+track+"&onlyIDs="+onlyIDs,
                    type: "GET",
                    datatype:"json",
                    async: false,    // default: true
                    success: function(result)
                    {
                        if(_DEBUG) console.log("annotationsGenes(): get annotations of genes done...");
                        //console.log(result);
                        response = result.response;
                    },
                    error: function()
                    {
                        if(_DEBUG) console.log("annotationsGenes(): annotationsGenes failed...");
                    }
                });
            return response;
        };


        Server.allAnnotationsGenes = function (callback, allPoints, types, window, align, onlyIDs,
                                               chromosomes, ws,
                                               annotations, numChromosome, startTime, numMatches)
        {
            // Initialize variables first
            if(typeof(annotations) === 'undefined')     annotations = "";
            if(typeof(numChromosome) === 'undefined')   numChromosome = 1;
            if(typeof(startTime) === 'undefined')       startTime = new Date();
            if(typeof(numMatches) === 'undefined')      numMatches = 1;


            // Show the image of "loading...".
            if(_loadingImage != null && annotations == "") { _loadingImage.css("top", 0).show(); }


            // Calculate track, points of this track and number of matches
            var track  = chromosomes[numChromosome-1];
            var points = JSON.parse(allPoints[track]);
            for( var j=0;j<points.length;j++) points[j]*=ws;
            numMatches += points.length;
            if(_DEBUG) console.log(points.length+" matches in track "+track);

            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/annotations?user="+_user+"&password="+_password+
                    "&positions=["+points+"]&types="+types+"&window="+window+"&align=\""+align+"\"&track="+track+"&onlyIDs="+onlyIDs,
                    type: "GET",
                    datatype: "json"
                });


            $.when(requestAJAX)
                .done(function(result)
                {
                    if(annotations != "" && result.response != "") annotations += ",";
                    annotations += result.response;
                    numChromosome++;

                    // While the number of chromosomes is less than or equal, we continue to get more annotations
                    if(numChromosome <= chromosomes.length)
                    {
                        Server.allAnnotationsGenes(getEnrichment, allPoints, "[\"gene\"]", window, "left", "True",
                            chromosomes, GVB_GLOBAL.ws,
                            annotations, numChromosome, (new Date()-startTime), numMatches);
                    }
                    else
                    {
                        // Hide the image of "loading...".
                        if(_loadingImage != null) { _loadingImage.css("display", "none"); }

                        if(_DEBUG) console.log("Total number of matches: "+numMatches);
                        if(_DEBUG) console.log("allAnnotationsGenes(): get all annotations done...");
                        if(_DEBUG) console.log("Time spent getting all annotations: "+ (new Date()-startTime)+"ms");

                        callback(annotations);
                    }
                })
                .fail(function()
                {
                    if(_DEBUG) console.log("allAnnotationsGenes(): get all annotations failed...");
                });
        };


        Server.enrichmentGO = function (callback, annotations, correction, alpha)
        {
            var startTime = new Date();

            // Show the image of "loading...".
            if(_loadingImage != null) { _loadingImage.css("top", 0).show(); }


            var arrayAnnotations = annotations;
            arrayAnnotations = arrayAnnotations.split(',');
            console.log("Number of gene annotations:"+arrayAnnotations.length);


            annotations = "[\""+annotations+"\"]";
            annotations = annotations.replace(/,/g,"\",\"");  //replaces ',' by '","'


            var requestAJAX = $.ajax(
                {
                    url: _serverPath+"/enrichmentGO?user="+_user+"&password="+_password+"&annotations="+annotations+"&correction="+correction+"&alpha="+alpha,
                    type: "GET",
                    datatype: "json"
                });

            $.when(requestAJAX)
                .done(function(result)
                {
                    // Convert to array...
                    var response = [];
                    response = result.response;

                    // Hide the image of "loading...".
                    if(_loadingImage != null) { _loadingImage.css("display", "none"); }

                    if(_DEBUG) console.log("enrichmentGO(): enrichmentGO done...");
                    if(_DEBUG) console.log("Time spent Enrichment analysis: "+ (new Date()-startTime)+"ms");

                    callback(response);

                })
                .fail(function()
                {
                    if(_DEBUG) console.log("enrichmentGO(): enrichmentGO failed...");
                });
        };


        Server.getPartSeq = function (track, startSeq, endSeq)
        {
            var response = [];
            $.ajax(
                {
                    url: _serverPath+"/getPartSeq?user="+_user+"&password="+_password+"&track="+track+"&start="+startSeq+"&end="+endSeq,
                    type: "GET",
                    datatype: "json",
                    async: false,    // default: true
                    success: function(result)
                    {
                        response.partSeq = result.partSeq;
                    },
                    error: function()
                    {
                        if(_DEBUG) console.log("getPartSeq(): getPartSeq failed...");
                    }
                });
            return response;
        };


        /**
         * Calls the REST service available to retrieve a sequence of nucleotides
         * @param start
         * @param end
         * @param track
         * @returns sequence for that interval
         */
        Server.nucleotides = function (start,end,track)
        {
            var response=[];
            $.ajax(
                {
                    url: _serverPath+"/nucleotides?user="+_user+"&password="+_password+"&start="+start+"&end="+end+"&track="+track,
                    type: "GET",
                    datatype:"json",
                    async: false,    // default: true
                    success: function(result)
                    {
                        if(_DEBUG) console.log("nucleotides(): done");
                        response.seq=result.response; // this is only a sample, as it is too large to show as a whole and to send via REST
                    },
                    error: function()
                    {
                        if(_DEBUG) console.log("preprocress(): discretization failed...");
                    }
                });
            return response;
        };


        /**
         * Calls the REST service available to retrieve a sequence of nucleotides
         * @param positions: array with numerical starting positions
         * @param size: length from the starting points
         * @param track: chromosome or track where the sequences must be taken from
         * @returns sequence for that interval
         */
        Server.nucProfile = function (positions,size,track)
        {
            var response=[];
            $.ajax(
                {
                    url: _serverPath+"/nucProfile?user="+_user+"&password="+_password+"&positions="+positions+"&size="+size+"&track="+track,
                    type: "GET",
                    datatype:"json",
                    async: false,    // default: true
                    success: function(result)
                    {
                        if(_DEBUG) console.log("nucProfile(): done");
                        response.seqs=result.seqs;
                        response.profile=result.profile;
                        response.consensus=result.consensus;
                        //The following three may be undefined if the number of positions is larger than 50-100
                        response.alignment=result.alignment;
                        response.aconsensus=result.aconsensus;
                        response.aprofile=result.aprofile;
                        //these ones are related to motif finding
                        response.motifs=result.motifs;
                        response.locations=result.locations;//motifs locations inside seqs
                        response.motifConsensus=result.motifConsensus;
                        response.motifProfile=result.motifProfile;
                    },
                    error: function()
                    {
                        if(_DEBUG) console.log("preprocress(): discretization failed...");
                    }
                });
            return response;
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
            msg   || ( msg = 'This is not an error. This is just to abort javascript' );

            // Hide the image of "loading...".
            if(_loadingImage != null) { _loadingImage.css("display", "none"); }

            alert(msg);

            throw "ERROR GBV: "+msg;
        }

        return Server;
    }

})(window);

