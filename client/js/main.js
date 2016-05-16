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
    filename: null,         // filename: name of file
    chromosomes: null,      // chromosomes: array with the names of all chromosomes
    track: null,            // track: name of selected chromosome
    ws: 30,                 // window size: discrete to real ratio
    nb: 5,                  // num bins
    maxSize: 400000         // maximum number of normalized data to store
};


// CHECK FILE (AND UPLOAD FILE)
////////////////////////////////////////////
function checkFile(file)
{
    destroyAll(false);

    // We ensure that there is a selected file
    if(typeof(file) != 'undefined')
    {
        if (DEBUG_GBV) console.log("\n----- UPLOAD FILE -----");
        if (DEBUG_GBV) console.log("Name of file: " + file.name);


        GVB_GLOBAL.filename = file.name;

        var forceReload = false;
        if ($("#reload").is(':checked'))
            forceReload = true;

        // We try to communicate with the server, upload  file (if necessary)
        Server.sendFile(preprocessing, file, forceReload);
    }
}


// PREPROCESSING: NORMALIZE, STATISTICAL DESCRIPTORS, DISCRETIZE...
/////////////////////////////////////////////////////////////////////////
function preprocessing(chromosome)
{
    if (DEBUG_GBV) console.log("\n----- PREPROCESSING -----");

    if(typeof(chromosome) === 'undefined')
    {
        GVB_GLOBAL.track = "None";
        if (DEBUG_GBV) console.log("chromosome: " + GVB_GLOBAL.track + " (first chromosome found)");
        Server.preprocess(drawingFirstDataLine, GVB_GLOBAL.filename, GVB_GLOBAL.track, GVB_GLOBAL.ws, GVB_GLOBAL.nb, GVB_GLOBAL.maxSize);
    }
    else
    {
        GVB_GLOBAL.track = chromosome;
        if (DEBUG_GBV) console.log("chromosome: " + GVB_GLOBAL.track);
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
    var chromosome          = GVB_GLOBAL.track;
    var numNucleotidesDraw  = globalDL1.cv.dim.width; // because the scale is 1:1

    dataLine_1_drawPoints(result.points, GVB_GLOBAL.chromosomes, chromosome, result.sizePattern, numNucleotidesDraw);


    // ENRICHMENT
    //----------------------------------
    if (DEBUG_GBV) console.log("\n----- ENRICHMENT -----");
    getAllAnnotations(result.points, result.sizePattern);
}


// GET ALL ANNOTATIONS
////////////////////////////////
function getAllAnnotations(allPoints, sizePattern)
{
    var chromosomes = GVB_GLOBAL.chromosomes;
    var ws          = GVB_GLOBAL.ws;

    Server.allAnnotationsGenes(getEnrichment, allPoints, "[\"gene\"]", sizePattern*ws, "left", "True",
                                chromosomes, ws);
}


// GET ENRICHMENT
////////////////////////////////
function getEnrichment(gis)
{
    /*
        // Jonatan: esta parte habría que implementarla en Server.allAnnotationsGenes, para que devolviera
        // también "annotations", propagándolo hacia dataLine_1_drawEnrichment().

        //Will solve posterior annotation queries, but it's very time consuming (increases from 1 to 10
        var gis=""//genes of interest
        var annotations={}//for each position, the genes in it
        for(var i=0;i<processedData.chromosomes.length;i++) {
            var points=JSON.parse(result.points[processedData.chromosomes[i]]);
            for( var j=0;j<points.length;j++)
                points[j]*=ws;
            console.log(points.length+" matches in track "+processedData.chromosomes[i]);
            numMatches+=points.length
            //annotations[processedData.chromosomes[i]]=Server.annotationsGenes("[" + points + "]", "[\"gene\"]", result.sizePattern*ws, "left", processedData.chromosomes[i], "False");
            annotations[processedData.chromosomes[i]]=Server.annotationsGenes("[" + points + "]", "[\"any\"]", globalDL2.dim.width, "center", processedData.chromosomes[i], "False");//this here might be too burdening
            gis+=Server.annotationsGenes("[" + points + "]", "[\"gene\"]", result.sizePattern*ws, "left", processedData.chromosomes[i], "True");
        }
        gis=gis.replace(/,$/, "");
        setAnnotations(gis,annotations)
    */

    Server.enrichmentGO(drawEnrichment, gis, "fdr", 0.01)
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
        $('#imagesChromosomes').append('<img style="margin-top:15px;margin-right:5px"'+
            'id="'+chromosomes[i]+'" class="image-chromosome" data-chromosome="'+chromosomes[i]+'" '+
            'src="images/chromosome.png" height="24px" width="24px">');
    }

    $(".image-chromosome").bind( "click", function()
    {
        var chromosomeName = $(this).data('chromosome');

        $(".image-chromosome").attr("src", "images/chromosome.png");
        $(this).attr("src", "images/chromosome_selected.png");

        preprocessing(chromosomeName);
    });

    // Tooltip to display information chromosome
    $(".image-chromosome").hover(
        // Move the mouse within the image.
        function()
        {
            var chromosomeName = $(this).data('chromosome');

            var pName = $('<p class="chromosome-name"></p>').text(chromosomeName);

            console.log("Chromosome name: ",chromosomeName);
            console.log(pName);

            $('<p class="chromosome-tooltip"></p>')
                .appendTo('body')
                .append(pName)
                .fadeIn('slow');
        },
        // Move the mouse away from the image
        function()
        {
            $('.chromosome-tooltip').remove();
        }
    ).mousemove(
        function(e)
        {
            var mouse_x0 = e.pageX + 5;
            var mouse_y0 = e.pageY + 5;
            $('.chromosome-tooltip').css({ top: mouse_y0, left: mouse_x0 })
        }
    );


    $("#"+chromosomes[0]).attr("src", "images/chromosome_selected.png");

}

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


    // Clear console and files configuration
    if(clear && DEBUG_GBV)
    {
        // Reset input files
        $("#files").val('');

        console.log(new Array(15).join("\n"));
        console.log("Reset all (with File APIs)...");
    }
}
