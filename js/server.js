


var serverPath="http://127.0.0.1:5000/";


// To try calls:
// curl -i -H "Accept: application/json" -H "Content-Typ: application/json" -X GET http://localhost:5000/testUpload?user=jpiri27&password=ninguna&filename=Mei3h_center_wl-peque2.wig&md5=51d4ad140a8346ceae6bca385058871b


/**
 * Uploads the file f to the server. Posterior calls to REST will be by the path of the file
 * @param file
 * @param hashMD5
 */
// NOTE: Remember REST Flask calls require enable CORS in the browser!!!

function sendFile(file, hashMD5)
{
    var response="";
    $.ajax(
    {
        url: serverPath+"testUpload?user="+user+"&password="+password+"&filename="+file.name+"&md5="+hashMD5,
        type: "GET",
        datatype:"json",
        async: false,
        success: function(result)
        {
            if(DEBUG) console.log("sendFile(): uploaded? "+result.response);
            response=result.response;
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
function preprocess(filename,track,ws,nb,maxSize)
{
    var response=[];
    $.ajax(
        {
            url: serverPath+"preprocess?user="+user+"&password="+password+"&filename="+filename+"&track="+track+"&windowSize="+ws+"&numBins="+nb+"&maxSize="+maxSize,
            type: "GET",
            datatype:"json",
            async: false,    // default: true
            success: function(result)
            {
                if(DEBUG) console.log("preprocress(): discretization done...");
                response.seq=result.seq; // this is only a sample, as it is too large to show as a whole and to send via REST
                response.max=result.maximum;
                response.min=result.minimum;
                response.mean=result.mean;
                response.stdev=result.stdev;
                response.dseq=result.dseq;
                response.fullLength=result.fullLength;
            },
            error: function(textStatus, errorThrown)
            {
                if(DEBUG) console.log("preprocress(): discretization failed...");
            }
        });
    return response;
}



function search(pattern,d)
{
    var response=[];
    $.ajax(
        {
            url: serverPath+"search?user="+user+"&password="+password+"&pattern="+pattern+"&d="+d,
            type: "GET",
            datatype:"json",
            async: false,    // default: true
            success: function(result)
            {
                if(result.response != "error")
                {
                    if(DEBUG) console.log("search(): search done...");

                    // Convert to array...
                    response.points = JSON.parse(result.response);
                }
                else
                {
                    console.log("search(): "+result.msg);
                }

            },
            error: function(textStatus, errorThrown)
            {
                console.log("search(): search failed...");
            }
        });
    return response;
}