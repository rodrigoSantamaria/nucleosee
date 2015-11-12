


var serverPath="http://127.0.0.1:5000/";


// To try calls:
// curl -i -H "Accept: application/json" -H "Content-Typ: application/json" -X GET http://localhost:5000/testUpload?filename=Mei3h_center_wl-peque2.wig


/**
 * Uploads the file f to the server. Posterior calls to REST will be by the path of the file
 * @param DEBUG
 * @param file
 * @param hashMD5
 */
// NOTE: Remember REST Flask calls require enable CORS in the browser!!!
function sendFile(DEBUG, user, password, filename, hashMD5)
{
    var ret="";
    $.ajax(
    {
        url: serverPath+"testUpload?user="+user+"&password="+password+"&filename="+filename+"&md5="+hashMD5,
        type: "GET",
        datatype:"json",
        async: false,
        success: function(result)
        {
            if(DEBUG) console.log("sendFile(): uploaded? "+result.response);
            ret=result.response;
        }
    });
    // only in this case we upload it
    if(ret=="outdated version" || ret=="not found")
    {
        // The FormData object lets you compile a set of key/value pairs to send using XMLHttpRequest (AJAX)
        var fd = new FormData();
        fd.append('file',file);

        $.ajax(
        {
            url: serverPath + "upload?user=jpiriz&password=ninguna",
            type: "POST",
            data: fd,
            processData: false, // tell jQuery not to process the data
            contentType: false, // tell jQuery not to set contentType
            async: false,
            success: function (result)
            {
                if(DEBUG) console.log("sendFile(): uploaded to "+result.path);
            }
        });
    }
}


/**
 * Calls the REST service available to a summary and characterization of data
 * @param DEBUG
 * @param path  local to the server where the data are present (usually a .wig file)
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
function preprocess(DEBUG,user,password,filename,track,ws,nb,maxSize)
{
    var result=[];
    $.ajax(
        {
            url: serverPath+"preprocess?user="+user+"&password="+password+"&filename="+filename+"&track="+track+"&windowSize="+ws+"&numBins="+nb+"&maxSize="+maxSize,
            type: "GET",
            datatype:"json",
            async: false,    // default: true
            success: function(response)
            {
                if(DEBUG) console.log("preprocress(): discretization done...");
                result.seq=response.seq; // this is only a sample, as it is too large to show as a whole and to send via REST
                result.max=response.maximum;
                result.min=response.minimum;
                result.mean=response.mean;
                result.stdev=response.stdev;
                result.dseq=response.dseq;
                result.fullLength=response.fullLength;
            },
            error: function(textStatus, errorThrown)
            {
                if(DEBUG) console.log("preprocress(): discretization failed...");
            }
        });
    return result;
}



function search(DEBUG,pattern,d)
{
    var result=[];
    $.ajax(
        {
            url: serverPath+"search?user=jpiriz&password=ninguna&pattern="+pattern+"&d="+d,
            type: "GET",
            datatype:"json",
            async: false,    // default: true
            success: function(response)
            {
                if(DEBUG) console.log("search(): search done...");

                console.log(response);

            },
            error: function(textStatus, errorThrown)
            {
                if(DEBUG) console.log("search(): search failed...");
            }
        });
    return result;
}