/*
  ┌────────────────────────────────────────────────────────────┐
  │ main.js                                                    │
  ├────────────────────────────────────────────────────────────┤
  │ Description:                                               │
  └────────────────────────────────────────────────────────────┘
*/


var user="jpiriz";
var password="ninguna";

var DEBUG_GBV = true;



//************************************************************************
//************************************************************************
// ---------------------------- Main code --------------------------------
//************************************************************************
//************************************************************************



function main(file, hashMD5)
{
    var track=0;         // track: number of cromosome
    var ws=30;           // window size: discrete to real ratio
    var nb=5;            // num bins
    var maxSize=400000;  // maximum number of normalized data to store


    destroyAll(false);

    Server.connect(DEBUG_GBV, user, password);

    var processedData=uploadFileAndPreprocess(file, hashMD5, track, ws, nb, maxSize);

    drawing(processedData, ws, maxSize);
}


// ANALYSIS: UPLOAD FILE AND PREPROCESSING
////////////////////////////////////////////
function uploadFileAndPreprocess(file, hashMD5, track, ws, nb, maxSize)
{
    // READ FILE
    //------------
    if(DEBUG_GBV) console.log("\n----- UPLOAD FILE -----");
    if(DEBUG_GBV) console.log("Name of file: "+file.name);

    var startTime=new Date();
    Server.sendFile(file, hashMD5);    // passing it to the server side (best solution for >1MB files)
    if (DEBUG_GBV) console.log("Time spent sending: " + (new Date() - startTime) + "ms");


    // PREPROCESSING: NORMALIZE, STATISTICAL DESCRIPTORS, DISCRETIZE
    //----------------------------------------------------------------
    if (DEBUG_GBV) console.log("\n----- PREPROCESSING -----");

    var fastMode=0;
    if($('#fastMode').is(":checked")) fastMode=1;

    startTime = new Date();
    var processedData = Server.preprocess(fastMode, file.name, track, ws, nb, maxSize);
    if (DEBUG_GBV) console.log("Time spent preprocessing: " + (new Date() - startTime) + "ms");

    return processedData;
}


// DRAWING
////////////////////////////////
function drawing(processedData, ws, maxSize)
{
    if(DEBUG_GBV) console.log("\n----- DRAWING -----");

    var seqServer=processedData.seq;    // just a sampling of about 100K of the original full length sequence
    var fullLength=processedData.fullLength;
    if(DEBUG_GBV) console.log("Length of seqServer:"+seqServer.length+" (full length seq="+fullLength+")");


    // DATALINE 1 (preprocessed data)
    //----------------------------------
    var startTime=new Date();
    dataLine1(seqServer,0,fullLength,fullLength, maxSize, processedData.mean, processedData.stdev, ws);
    if(DEBUG_GBV) console.log("Time spent dataLine1: "+ (new Date()-startTime)+"ms");

}


// SEARCH POINTS
////////////////////////////////
function searchPoints()
{
    if(globalDL1.drawn)
    {
        var pattern = $('#patternSearch').val();
        var d       = $('#dSearch').val();

        // DATALINE 2: DRAW POINTS
        //----------------------------------
        var numNucleotidesDraw = globalDL1.width;

        var startTime=new Date();
        var result=Server.search(pattern,d);
        result.points = [6033]; // TODO: quitar esta línea (points)
        dataLine1_drawPoints(result.points, result.sizePattern, numNucleotidesDraw);
        if(DEBUG_GBV) console.log("Time spent search: "+ (new Date()-startTime)+"ms");


        // DATALINE 3: GET ANNOTATIONS
        //------------------------------
        globalDL1.seqPoints = [180990]; // TODO: quitar esta línea (points)
        Server.annotationsGenes(globalDL3, "["+globalDL1.seqPoints+"]","[\"any\"]",globalDL1.width*2);
    }
}




// DESTROY ALL
////////////////////////////////
function destroyAll(clear)
{
    // Reset input files
    $("#files").val('');

    // Empty all SVG images
    var images = ["#lineSeq", "#lineSeq2"];

    for(var i=0;i<images.length;i++)
    {
        if($(images[i]).html() != "")
            $(images[i]).empty();
    }

    // Clear console
    if(clear && DEBUG_GBV)
    {
        console.log(new Array(15).join("\n"));
        console.log("Reset all (with File APIs)...");
    }
}


