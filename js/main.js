



//************************************************************************
//************************************************************************
// ---------------------------- Main code --------------------------------
//************************************************************************
//************************************************************************



function main(user, password, file, hashMD5)
{
    destroyAll(false);

    var desc=uploadFileAndProcessing(true, user, password, file, hashMD5);

    drawing(true, desc);
}


// ANALYSIS: UPLOAD FILE AND PROCESSING
/////////////////////////////////////////
function uploadFileAndProcessing(DEBUG, user, password, file, hashMD5)
{
    // READ FILE
    //------------
    if(DEBUG) console.log("\n----- UPLOAD FILE -----");
    if(DEBUG) console.log("Name of file: "+file.name);

    var startTime=new Date();
    sendFile(DEBUG, user, password, file, hashMD5);    // passing it to the server side (best solution for >1MB files)
    if(DEBUG) console.log("Time spent sending: "+ (new Date()-startTime)+"ms");


    // PREPROCESSING: NORMALIZE, STATISTICAL DESCRIPTORS, DISCRETIZE
    //----------------------------------------------------------------
    if(DEBUG) console.log("\n----- PROCESSING -----");

    var track=0;
    var ws=150;         // window size: discrete to real ratio
    var nb=5;           // num bins
    var maxSize=100000;  // maximum number of normalized data to store

    startTime=new Date();
    var desc=preprocess(DEBUG, user, password, file.name, track, ws, nb, maxSize);
    if(DEBUG) console.log("Time spent preprocessing: "+ (new Date()-startTime)+"ms");

    return desc;
}




// DRAWING
////////////////////////////////
function drawing(DEBUG, desc)
{
    if(DEBUG) console.log("\n----- DRAWING -----");

    var seq=desc.seq;    // just a sampling of about 100K of the original full length sequence
    var min=desc.min;
    var max=desc.max;
    var mean=desc.mean;
    var stdev=desc.stdev;
    var disc=desc.dseq; // whole seq compression
    var fullLength=desc.fullLength;
    if(DEBUG) console.log("Length of seq:"+seq.length+" (full length="+fullLength+")");


    // Dimensions of screen
    var graphHeight=200;
    var graphWidth=screen.width;
    var startX=50;
    var startY=50;

    // DATA LINE (preprocessed data)
    //--------------------------------
    var startTime=new Date();
    var rDataLine = dataLine(DEBUG, seq, 0, seq.length, mean, stdev, graphHeight, graphWidth, startX, startY);
    if(DEBUG) console.log("Time spent dataline: "+ (new Date()-startTime)+"ms");








    var dataPoints=[];
    dataPoints.push({pos: (rDataLine.data)[900].pos, value: (rDataLine.data)[900].value});
    dataPoints.push({pos: (rDataLine.data)[600].pos, value: (rDataLine.data)[600].value});
    dataPoints.push({pos: (rDataLine.data)[400].pos, value: (rDataLine.data)[400].value});
    dataPoints.push({pos: (rDataLine.data)[300].pos, value: (rDataLine.data)[300].value});

    drawPoints(true, rDataLine.svg,dataPoints,rDataLine.x,rDataLine.y,rDataLine.window,
                seq, mean, stdev, graphHeight, graphWidth, startX, startY);

}


// DESTROY ALL
////////////////////////////////
function destroyAll(clear)
{
    // Reset input files
    $("#files").val('');

    // Empty all SVG images
    if($('#lineSeq').html() != "")
        $("#lineSeq").empty();

    if($('#lineSeq2').html() != "")
        $("#lineSeq2").empty();

    if(clear)
    {
        // Clear console
        console.log(new Array(15).join("\n"));
        console.log("Reset all (with File APIs)...");
    }
}


function seachPoints()
{

    var pattern = $('#patternSearch').val();
    var d       = $('#dSearch').val();

    search(true,pattern,d);

}
