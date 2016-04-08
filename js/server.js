/*
 ┌────────────────────────────────────────────────────────────┐
 │ server.js                                                  │
 ├────────────────────────────────────────────────────────────┤
 │ Description:                                               │
 └────────────────────────────────────────────────────────────┘
 */


// To try calls:
// curl -i -H "Accept: application/json" -H "Content-Typ: application/json" -X GET http://localhost:5000/testUpload?user=jpiri27&password=ninguna&filename=Mei3h_center_wl-peque2.wig&md5=51d4ad140a8346ceae6bca385058871b


(function(window){

    //I recommend this
    'use strict';

    //define globally if it doesn't already exist
    if(typeof(Server) === 'undefined') {
        window.Server = defineServerLibrary();
    }
    else {
        console.log("Library already defined.");
    }



    function defineServerLibrary()
    {
        var Server = {};

        var _serverPath = "";
        var _user       = "";
        var _password   = "";
        var _DEBUG      = false;


        Server.connect = function (DEBUG, user, passwordserver, serverPath)
        {
            serverPath   || ( serverPath = "http://127.0.0.1:5000/" );

            _serverPath = serverPath;
            _user       = user;
            _password   = password;
        };


        /**
         * Uploads the file f to the server. Posterior calls to REST will be by the path of the file
         * @param file
         * @param forceReload
         */
        // NOTE: Remember REST Flask calls require enable CORS in the browser!!!
        Server.sendFile = function (file, forceReload)
        {
            var response="";
            $.ajax(
                {
                    //url: _serverPath+"testUpload?user="+_user+"&password="+_password+"&filename="+file.name+"&md5="+forceReload,
                    url: _serverPath+"testUpload?user="+_user+"&password="+_password+"&filename="+file.name+"&forceReload="+forceReload,
                    type: "GET",
                    datatype:"json",
                    async: false,
                    success: function(result)
                    {
                        if(_DEBUG) console.log("sendFile(): uploaded? "+result.response);
                        response=result.response;
                    },
                    error: function(XMLHttpRequest, textStatus, errorThrown)
                    {
                        if(_DEBUG) console.log("sendFile(): testUpload failed...");
                        javascript_abort(XMLHttpRequest.statusText);
                    }


                });
            // only in this case we upload it
            if(response=="outdated version" || response=="not found")
            {
                // The FormData object lets you compile a set of key/value pairs to send using XMLHttpRequest (AJAX)
                var fd = new FormData();
                fd.append('file',file);

                $.ajax(
                    {
                        url: _serverPath+"upload?user="+_user+"&password="+_password,
                        type: "POST",
                        data: fd,
                        processData: false, // tell jQuery not to process the data
                        contentType: false, // tell jQuery not to set contentType
                        async: false,
                        success: function (result)
                        {
                            if(_DEBUG) console.log("sendFile(): uploaded to "+result.path);
                        },
                        error: function()
                        {
                            if(_DEBUG) console.log("sendFile(): upload failed...");
                        }
                    });
            }
        };


        /**
         * Calls the REST service available to a summary and characterization of data
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
        Server.preprocess = function (filename,track,ws,nb,maxSize)
        {
            var response=[];
            $.ajax(
                {
                    url: _serverPath+"preprocess?user="+_user+"&password="+_password+"&filename="+filename+"&track="+track+"&windowSize="+ws+"&numBins="+nb+"&maxSize="+maxSize,
                    type: "GET",
                    datatype:"json",
                    async: false,    // default: true
                    success: function(result)
                    {
                        if(_DEBUG) console.log("preprocress(): discretization done...");
                        response.seq=result.seq; // this is only a sample, as it is too large to show as a whole and to send via REST
                        response.max=result.maximum;
                        response.min=result.minimum;
                        response.mean=result.mean;
                        response.stdev=result.stdev;
                        response.dseq=result.dseq;
                        response.fullLength=result.fullLength;
                        response.chromosomes=result.chromosomes;
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
                    url: _serverPath+"nucleotides?user="+_user+"&password="+_password+"&start="+start+"&end="+end+"&track="+track,
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
                    url: _serverPath+"nucProfile?user="+_user+"&password="+_password+"&positions="+positions+"&size="+size+"&track="+track,
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
                    },
                    error: function()
                    {
                        if(_DEBUG) console.log("preprocress(): discretization failed...");
                    }
                });
            return response;
        };

        Server.search = function (pattern,d)
        {
            var response=[];
            pattern=pattern.replace(/\+/g, "%2B");
            console.log("pattern is "+pattern);
            $.ajax(
                {
                    url: _serverPath+"search?user="+_user+"&password="+_password+"&pattern="+pattern+"&d="+d,
                    type: "GET",
                    datatype:"json",
                    async: false,    // default: true
                    success: function(result)
                    {
                        if(result.response != "error")
                        {
                            if(_DEBUG) console.log("search(): search done...");

                            // Convert to array...
                            //response.points = JSON.parse(result.points);
                            response.points=result.points;
                            response.sizePattern = result.sizePattern;
                        }
                        else
                        {
                            if(_DEBUG) console.log("search(): "+result.msg);
                        }
                    },
                    error: function()
                    {
                        if(_DEBUG) console.log("search(): search failed...");
                    }
                });
            return response;
        };


        Server.getPartSeq = function (startSeq,endSeq, track)
        {
            var response=[];
            $.ajax(
                {
                    url: _serverPath+"getPartSeq?user="+_user+"&password="+_password+"&start="+startSeq+"&end="+endSeq+"&track="+track,
                    type: "GET",
                    datatype:"json",
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


        Server.annotationsGenes = function (points,types,window, align, track, onlyIDs)
        {
            var response=[];
            $.ajax(
                {
                    url: _serverPath+"annotations?user="+_user+"&password="+_password+"&positions="+points+"&types="+types+"&window="+window+"&align=\""+align+"\"&track="+track+"&onlyIDs="+onlyIDs,
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

        Server.enrichment = function (annotations, correction, alpha)
        {
            annotations="[\""+annotations+"\"]";
            annotations=annotations.replace(/,/g,"\",\"");
            console.log("Number of gene annotations:"+annotations.length)
            var response=[];
            $.ajax(
                {
                    url: _serverPath+"enrichmentGO?user="+_user+"&password="+_password+"&annotations="+annotations+"&correction="+correction+"&alpha="+alpha,
                    type: "GET",
                    datatype:"json",
                    async: false,    // default: true
                    success: function(result)
                    {
                        if(_DEBUG) console.log("enrichmentGO(): done");
                        response = result.response;

                    },
                    error: function()
                    {
                        if(_DEBUG) console.log("annotationsGenes(): annotationsGenes failed...");
                    }
                });
            return response;
        };

        function javascript_abort(msg)
        {
            msg   || ( msg = 'This is not an error. This is just to abort javascript' );

            throw "ERROR GBV: "+msg;
        }

        return Server;
    }

})(window);


































function SERVER_getPartSeq(startSeq,endSeq)
{
    var response=[];
    $.ajax(
        {
            url: serverPath+"getPartSeq?user="+user+"&password="+password+"&start="+startSeq+"&end="+endSeq,
            type: "GET",
            datatype:"json",
            async: false,    // default: true
            success: function(result)
            {
                response.partSeq = result.partSeq;
            },
            error: function()
            {
                console.log("getPartSeq(): getPartSeq failed...");
            }
        });
    return response;
}


function SERVER_annotationsGenes(points,types,window)
{
    var response=[];
    $.ajax(
        {
            url: serverPath+"annotations?user="+user+"&password="+password+"&positions="+points+"&types="+types+"&window="+window,
            type: "GET",
            datatype:"json",
            async: true,    // default: true
            success: function(result)
            {
                console.log(result);
            },
            error: function()
            {
                console.log("annotationsGenes(): annotationsGenes failed...");
            }
        });
    return response;
}