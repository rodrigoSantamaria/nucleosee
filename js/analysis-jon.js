
/** TODO: all this should go to python on the server side -> undergoing work via analysis.py methods
 * Created by rodri on 27/01/14.
 */
var serverPath="http://127.0.0.1:5000/" //it will be changed to the IP where the service will be hosted
var filePath=""

/**
 * Uploads the file f to the server. Posterior calls to REST will be by the path of the file
 * @param f
 * @returns {string} path to the file in the server
 */
//NOTE: Remember REST Flask calls require enable CORS in the browser!!! <----
function sendFile(DEBUG, f)
    {
    var fd = new FormData();
    fd.append( 'file', f );
    //TODO: before posting, which is costly if the file is large, we can check if it's already uploaded?
        //  At least for tests, this could prove useful
    var ret="";
    // curl -i -H "Accept: application/json" -H "Content-Typ: application/json" -X GET http://localhost:5000/testUpload?filename=Mei3h_center_wl-peque2.wig
    $.ajax({
        url: serverPath+"testUpload?filename="+ f.name,
        type: "GET",
        datatype:"json",
        success: function(response)
        {
            if(DEBUG) console.log("sendFile(): uploaded? "+response.response);
            ret=response.response;
        },
        async: false
    });
    if(ret=="not found") { //only in this case we upload it
        $.ajax({
            url: serverPath + "upload",
            type: "POST",
            data: fd,
            processData: false,
            contentType: false,
            success: function (response)
            {
                if(DEBUG) console.log("sendFile(): uploaded to " + response.path);
                filePath = response.path;
                ret = response.path;
            },
            async: false
        });
        return ret;
        }
    else
        return ret;
    }


/**
 * Calls the REST service available to a summary and characterization of data
 * @param path  local to the server where the data are present (usually a .wig file)
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
function preprocess(DEBUG, path,ws,nb,maxSize)
    {//TODO: save after preprocessing in the server and re-preprocess if not
    var na=[];
    $.ajax({
        url: serverPath+"preprocess?path="+path+"&windowSize="+ws+"&numBins="+nb+"&maxSize="+maxSize,
        type: "get",
        datatype:"json",
        success: function(response)
            {
            if(DEBUG) console.log("preprocress(): discretization done...");
            na.seq=response.result; //this is only a sample, as it is too large to show as a whole and to send via REST
            na.max=response.maximum;
            na.min=response.minimum;
            na.mean=response.mean;
            na.stdev=response.sdev;
            na.dseq=response.dseq;
            na.fullLength=response.fullLength;
            },
        async: false
        });
    return na;
    }



