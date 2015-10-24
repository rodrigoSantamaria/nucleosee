

// Dimensions of screen
screenWidth=screen.width;
startX=50;
startY=50;



//************************************************************************
//************************************************************************
//************************************************************************
// ---------------------------- Main code --------------------------------
//************************************************************************
//************************************************************************
//************************************************************************



// Check if there's support for file reading
if (window.File && window.FileReader && window.FileList && window.Blob)
    {
    console.log("Start all (with File APIs)...");
    this.handleEvent=function(evt)
    {
        switch(event.type)
            {
            case 'change':

                if ( $("#dna").length)   { $("#dna").remove(); var $div = $('<div />').appendTo('body');  $div.attr('id', 'dna');}
                if ( $("#dna2").length)  { $("#dna2").remove(); var $div = $('<div />').appendTo('body'); $div.attr('id', 'dna2');}

                var files = evt.target.files; // FileList object
                file=files[0]  //By now, just one file
                desc=uploadFileAndProcessing(true, file);

                drawing(true, desc);
            }
    }
    // Register the custom listener
    document.getElementById('files').addEventListener('change', this, false); //TODO: changing to a post method to upload files to server
    }
else
    {
    alert('The File APIs are not fully supported by your browser.');
    }


// UPLOAD FILE AND PROCESSING
////////////////////////////////
function uploadFileAndProcessing(DEBUG, file)
{
    // READ FILE
    //------------
    if(DEBUG) console.log("\n----- UPLOAD FILE -----");
    if(DEBUG) console.log("Name of file: "+file.name);

    startTime=new Date()
    var serverFilePath=sendFile(DEBUG, file);    // passing it to the server side (best solution for >1MB files)
    if(DEBUG) console.log("Path in remote: "+serverFilePath);
    if(DEBUG) console.log("Time spent sending: "+ (new Date()-startTime)+"ms");


    // PREPROCESSING: NORMALIZE, STATISTICAL DESCRIPTORS, DISCRETIZE
    //--------------------------------------------------------------
    if(DEBUG) console.log("\n----- PROCESSING -----");

    var ws=150;         // window size: discrete to real ratio
    var nb=5;           // num bins
    var maxSize=100000  // maximum number of normalized data to store

    startTime=new Date()
    var desc=preprocess(DEBUG, serverFilePath, ws, nb, maxSize);    // normalize+stats+discretize+suffix array build (next)
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
    var disc=desc.dseq; // whole seq compression
    var fullLength=desc.fullLength;

    if(DEBUG) console.log("Length of na:"+seq.length+" (full length="+fullLength+")");


    console.log(desc);

    graphHeight=100;
    graphWidth=screenWidth;

    //DATA LINE (preprocessed data)
    dataLine(seq, 0, seq.length, desc.mean, desc.stdev, graphHeight, graphWidth, startX, startY);


}

