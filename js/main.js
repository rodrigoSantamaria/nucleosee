/*
  ┌────────────────────────────────────────────────────────────┐
  │ main.js                                                    │
  ├────────────────────────────────────────────────────────────┤
  │ Description:                                               │
  └────────────────────────────────────────────────────────────┘
*/



var DEBUG_GBV = true;


var GVB_GLOBAL =
{
    filename: null,
    ws: 30,                 // window size: discrete to real ratio
    track: null,            // track: number of chromosome
    nb: 5,                  // num bins
    maxSize: 400000,        // maximum number of normalized data to store
    chromosomes: null
};


// CHECK FILE (AND UPLOAD FILE)
////////////////////////////////////////////
function checkFile(file)
{
    destroyAll(false);

    var forceReload = false;
    if($("#reload").is(':checked') )
        forceReload = true;


    if(DEBUG_GBV) console.log("\n----- UPLOAD FILE -----");
    if(DEBUG_GBV) console.log("Name of file: "+file.name);

    GVB_GLOBAL.filename = file.name;

    // We try to communicate with the server, upload  file (if necessary)
    Server.sendFile(preprocessing, file, forceReload);
}


// PREPROCESSING: NORMALIZE, STATISTICAL DESCRIPTORS, DISCRETIZE...
/////////////////////////////////////////////////////////////////////////
function preprocessing(chromosome)
{
    if (DEBUG_GBV) console.log("\n----- PREPROCESSING -----");

    if(typeof(chromosome) === 'undefined')
    {
        GVB_GLOBAL.track = "None";
        if (DEBUG_GBV) console.log("chromosome: "+GVB_GLOBAL.track+" (first chromosome found)");
        Server.preprocess(drawingFirstDataLine, GVB_GLOBAL.filename, GVB_GLOBAL.track, GVB_GLOBAL.ws, GVB_GLOBAL.nb, GVB_GLOBAL.maxSize);
    }
    else
    {
        GVB_GLOBAL.track = chromosome;
        if (DEBUG_GBV) console.log("chromosome: "+GVB_GLOBAL.track);
        Server.getTrack(drawingFirstDataLine, GVB_GLOBAL.track);
    }
}



// DRAWING: DATALINE 1
////////////////////////////////
function drawingFirstDataLine(processedData, chromosome)
{
    // We create icons chromosomes and bind the click event
    GVB_GLOBAL.chromosomes = processedData.chromosomes;
    if(chromosome == "None")
    {
        createIconsChromosomes(GVB_GLOBAL.chromosomes);
        GVB_GLOBAL.track    = processedData.chromosomes[0];
    }
    else
    {
        GVB_GLOBAL.track    = chromosome;
    }


    var seqServer       = processedData.seq;    // just a sampling of about 400K of the original full length sequence
    var fullLength      = processedData.fullLength;
    var mean            = processedData.mean;
    var stdev           = processedData.stdev;


    // DATALINE 1 (preprocessed data)
    //----------------------------------
    if(DEBUG_GBV) console.log("\n----- DRAWING: DATALINE 1 ("+GVB_GLOBAL.track+")-----");
    if(DEBUG_GBV) console.log("Length of seqServer:"+seqServer.length+" (full length seq="+fullLength+")");
    dataLine_1(GVB_GLOBAL.track, fullLength, seqServer, 0, fullLength, GVB_GLOBAL.maxSize, mean, stdev, GVB_GLOBAL.ws);
}





// SEARCH POINTS
////////////////////////////////
function searchPattern()
{
    if(globalDL1.drawn)
    {
        var pattern     = $('#patternSearch').val();
        var d           = $('#dSearch').val();

        if (DEBUG_GBV) console.log("\n----- SEARCH -----");
        Server.search(drawPoints, pattern,d);
    }
}







// DRAW POINTS ON DATALINE 1
////////////////////////////////
function drawPoints(result)
{
    // DATALINE 1: DRAW POINTS
    //----------------------------------
    var chromosome  = GVB_GLOBAL.track;
    var numNucleotidesDraw  = globalDL1.cv.dim.width; // because the scale is 1:1
    var points              = JSON.parse(result.points[chromosome]);

    dataLine_1_drawPoints(points, result.sizePattern, numNucleotidesDraw);



    // ENRICHMENT
    //----------------------------------
    if (DEBUG_GBV) console.log("\n----- ENRICHMENT -----");
    getAllAnnotations(result.points, result.sizePattern);




     //GET SEQUENCES, MOTIFS, (ALIGNMENT), CONSENSUS
     //NOTE: alignment takes more than 1s if there's >50 sequences! (using the fastest method: kalign)
     var start = new Date().getTime();
     for( var j=0;j<points.length;j++)
         points[j]*=GVB_GLOBAL.ws;
     var response=Server.nucProfile("["+points+"]",result.sizePattern*GVB_GLOBAL.ws,globalSeq.track); //TODO: by now only on this chromosome
     var end = new Date().getTime();
     console.log("Sequence analysis took: "+(end-start))
     setSequences(response)
}

// GET ALL ANNOTATIONS
////////////////////////////////
function getAllAnnotations(allPoints, sizePattern)
{
    var chromosomes = GVB_GLOBAL.chromosomes;

    Server.allAnnotationsGenes(getEnrichment, allPoints, "[\"gene\"]", sizePattern*GVB_GLOBAL.ws, "left", "True",
                                chromosomes, GVB_GLOBAL.ws);
}


// GET ENRICHMENT
////////////////////////////////
function getEnrichment(annotations)
{
    Server.enrichmentGO(drawEnrichment, annotations, "fdr", 0.01)
}


// DRAW ENRICHMENT
////////////////////////////////
function drawEnrichment(enrichment)
{
    dataLine_1_drawEnrichment(enrichment);
}







function createIconsChromosomes(chromosomes)
{
    $('#imagesChromosomes').empty();
    for(var i=0; i<chromosomes.length; i++)
    {
        $('#imagesChromosomes').append('<img style="margin-top:15px;margin-right:5px" id="'+chromosomes[i]+'" class="image-chromosome" data-chromosome="'+chromosomes[i]+'" src="images/chromosome.png" height="24px" width="24px">');
    }

    $(".image-chromosome").bind( "click", function()
    {
        var chromosomeImage = $(this).data('chromosome');

        $(".image-chromosome").css("background-color", "");
        $(this).css("background-color", "silver");

        preprocessing(chromosomeImage);
    });

    $("#"+chromosomes[0]).css("background-color", "silver");
}


// DESTROY ALL
////////////////////////////////
function destroyAll(clear)
{
    // Empty all SVG images (of the array)
    var images = ["lineSeq", "lineSeq2", "lineSeq3"];
    for(var i=0;i<images.length;i++)
    {
        var image = $("#"+images[i]);
        if(image.html() != "")
            image.empty();
    }

    // Empty all icons of chromosomes
    $("#imagesChromosomes").empty();

    // Clear console
    if(clear && DEBUG_GBV)
    {
        // Reset input files
        $("#files").val('');

        console.log(new Array(15).join("\n"));
        console.log("Reset all (with File APIs)...");
    }
}
