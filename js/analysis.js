
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
function sendFile(f)
    {
    var fd = new FormData();
    fd.append( 'file', f );
    //TODO: before posting, which is costly if the file is large, we can check if it's already uploaded?
        //  At least for tests, this could prove useful
    var ret="";
    $.ajax({
        url: serverPath+"testUpload?filename="+ f.name,
        type: "get",
        datatype:"json",
        success: function(response){
            console.log("file uploaded? "+response.response);
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
            success: function (response) {
                console.log("file uploaded to " + response.path);
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
 * Utility methods, deprecated by the use of REST services in python (see 'preprocess')
 * @param array
 */
function mean(v)
    {
    var sum=0
    for(i in v)
        sum+=parseFloat(v[i]);
    return sum/v.length;
    }

function sdev(v)
    {
    var m=mean(v);
    var sum=0;
    for(i in v)
        sum+=Math.pow(parseFloat(v[i])-m,2);
    return Math.sqrt(sum);
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
function preprocess(path,ws,nb,maxSize)
    {//TODO: save after preprocessing in the server and re-preprocess if not
    var na=[];
    $.ajax({
        url: serverPath+"preprocess?path="+path+"&windowSize="+ws+"&numBins="+nb+"&maxSize="+maxSize,
        type: "get",
        datatype:"json",
        success: function(response)
            {
            console.log("discretization done");
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

function normalizeZ(path)
    {
    var na=[];
    $.ajax({
        url: serverPath+"normalize?path="+path,
        type: "get",
        datatype:"json",
        success: function(response){
            console.log("discretization done");
            na=response.result;
            na.max=response.maximum;
            na.min=response.minimum;
            na.mean=response.mean;
            na.stdev=response.sdev;
        },
        async: false
    });
    return na;
    }

function mode(w)
    {
    var freqs={};
    console.log("this is : "+freqs[23]);
    for(i in w)
        {
        if(freqs[w[i]]==undefined)
            freqs[w[i]]=1;
        else
            freqs[w[i]]++;
        }
    var max=0;
    var maxi=0;
    for(i in freqs)
        if(freqs[i]>max)
            {
            max=freqs[i];
            maxi=i;
            }
    console.log("Most frequent element is "+maxi+", and appears "+freqs[maxi]+" times");
    return maxi;
    }

//As of Pearson's simplificaton
function skew(w)
    {
    //as in R package moments:
    var m=mean(w);
    var sumn=0;
    var sumd=0;
    for(i in w)
        {
        sumn+=Math.pow(parseFloat(w[i])-m,3);
        sumd+=Math.pow(parseFloat(w[i])-m,2);
        }
    sumn/= w.length;
    sumd/= w.length;
    console.log(sumn+", "+sumd);
    return sumn/Math.pow(sumd,3.0/2.0);
    }

function maximum(w)
{
    var max=w[0];
    for(i in w)
        if(max<w[i])    max=w[i];
    return max;
}
function minimum(w)
    {
    var min=w[0];
    for(i in w)
        if(min>w[i])    min=w[i];
    return min;
    }

function distribution(w, max, min, numBins)
    {
    var bin=(max-min)/numBins;
    var freqs=new Array(numBins);
    for(var i=0;i<numBins;i++)
        freqs[i]=0;
    for(i in w)
        {
        var index=Math.round((w[i]+Math.abs(min))/bin);
        freqs[index]+=1;
        }
    for(i in freqs)
        {
        freqs[i]/= w.length;
        }
    return freqs;
    }

/**
 * Given a sequence and a window size, it generates an alphabet a-e where a means that the average window value is
 * below mean-2*sdev, b < mean-sd, c<mean, d<mean+sd and e>mean+2*sd
 * @param seq
 * @param windowSize
 * @param mean
 * @param sdev
 */
function discretize(seq, windowSize, numBins)
    {
     var disc=[];
     $.ajax({
         url: serverPath+"discretize?numBins="+numBins+"&windowSize="+windowSize+"&seq="+seq,
         type: "get",
         datatype:"json",
         success: function(response){
         console.log("discretization done");
         disc=response.result;
         },
         async: false
         });
     return disc;
    }

function kurtosis(w)
    {
    var m=mean(w);
    var sum=0;

    //as implemented in R package 'moments'
    var sums=0;
    for(i in w)
        {
        sum+=Math.pow(parseFloat(w[i])-m,4);
        sums+=Math.pow(parseFloat(w[i])-m,2);
        }
    return w.length*sum/(sums*sums);
    }

// Returns a statistical description of the distribution at the window w
function statisticalDescription(w)
    {
    var description={}

    $.ajax({
        url: serverPath+"stats?seq="+w,
        type: "get",
        datatype:"json",
        success: function(response){
            console.log("stats computed");
            description.mean=response.mean;
            description.stdev=response.sdev;
            description.skew=0;//TODO
            description.k=0;//TODO
        },
        async: false
    });
//    description.mean=mean(w);
//    description.stdev=sdev(w);
//    description.k=kurtosis(w);
//    description.skew=skew(w);
    return description;
    }
