/*
  ┌────────────────────────────────────────────────────────────┐
  │ Main.js                                                    │
  ├────────────────────────────────────────────────────────────┤
  │ Description:                                               │
  └────────────────────────────────────────────────────────────┘
*/


var user="jpiriz";
var password="ninguna";

var DEBUG = true;



//************************************************************************
//************************************************************************
// ---------------------------- Main code --------------------------------
//************************************************************************
//************************************************************************



function main(file, hashMD5)
{
    var track=0;
    var ws=30;           // window size: discrete to real ratio
    var nb=5;            // num bins
    var maxSize=200000;  // maximum number of normalized data to store


    destroyAll(false);

    var desc=uploadFileAndProcessing(file, hashMD5, track, ws, nb, maxSize);

    drawing(desc, ws);
}


// ANALYSIS: UPLOAD FILE AND PROCESSING
/////////////////////////////////////////
function uploadFileAndProcessing(file, hashMD5, track, ws, nb, maxSize)
{
    // READ FILE
    //------------
    if(DEBUG) console.log("\n----- UPLOAD FILE -----");
    if(DEBUG) console.log("Name of file: "+file.name);

    var startTime=new Date();
    SERVER_sendFile(file, hashMD5);    // passing it to the server side (best solution for >1MB files)
    if(DEBUG) console.log("Time spent sending: "+ (new Date()-startTime)+"ms");


    // PREPROCESSING: NORMALIZE, STATISTICAL DESCRIPTORS, DISCRETIZE
    //----------------------------------------------------------------
    if(DEBUG) console.log("\n----- PREPROCESSING -----");

    startTime=new Date();
    var desc=SERVER_preprocess(file.name, track, ws, nb, maxSize);
    if(DEBUG) console.log("Time spent preprocessing: "+ (new Date()-startTime)+"ms");

    return desc;
}


// DRAWING
////////////////////////////////
function drawing(desc, ws)
{
    if(DEBUG) console.log("\n----- DRAWING -----");

    var seq=desc.seq;    // just a sampling of about 100K of the original full length sequence
    var fullLength=desc.fullLength;
    if(DEBUG) console.log("Length of seq:"+seq.length+" (full length="+fullLength+")");


    // DATA LINE (preprocessed data)
    //--------------------------------
    var startTime=new Date();
    dataLine1(seq, 100000, 0, fullLength, desc.mean, desc.stdev, ws);
    if(DEBUG) console.log("Time spent dataLine1: "+ (new Date()-startTime)+"ms");
}


// SEARCH POINTS
////////////////////////////////
function searchPoints()
{
    var pattern = $('#patternSearch').val();
    var d       = $('#dSearch').val();

    var startTime=new Date();
    var result = SERVER_search(pattern,d);
    dataLine1_drawPoints(result.points);
    if(DEBUG) console.log("Time spent search: "+ (new Date()-startTime)+"ms");
    //SERVER_annotationsGenes(globalDL1.seqPoints,["any"],globalDL1.width);
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
    if(clear && DEBUG)
    {
        console.log(new Array(15).join("\n"));
        console.log("Reset all (with File APIs)...");
    }
}


