/*
  ┌────────────────────────────────────────────────────────────┐
  │ main.js                                                    │
  ├────────────────────────────────────────────────────────────┤
  │ Description:                                               │
  └────────────────────────────────────────────────────────────┘
*/



var DEBUG_GBV = true;

var globalTime;

var GVB_GLOBAL =
{
    filename: null,         // filename: name of file
    chromosomes: null,      // chromosomes: array with the names of all chromosomes
    track: null,            // track: name of selected chromosome
    ws: 30,                 // window size: discrete to real ratio
    nb: 5,                  // num bins
    maxSize: 100000,        // maximum number of normalized data to store
    intersect: "soft",   // If true Hard instersects are considered (whole inclusion of pattern in genomic annotation for enrichment, etc.)
    geo: false,            // True if some genomic information (genes, UTRs, etc.) is used to filter out pattern matches
    softMutations: true,   // if true soft mutations (only 1-distance switch) is allowed to the patterns. E.g. "b" may change to "c" or "a" but not to "d" or "e".
    grid : false           // if true, a grid to show percentile regions is shown on both lane 1 and 2
};


// CHECK FILE (AND UPLOAD FILE)
////////////////////////////////////////////

function checkFile(file)
{
    /**
     * @typedef {Object} file
     * @property name
     */

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
    /**
     * @typedef {Object} processedData
     * @property seq
     * @property mean
     * @property stdev
     * @property dseq
     * @property fullLength
     * @property chromosomes
     */

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

    dataLine_1(GVB_GLOBAL.chromosomes, GVB_GLOBAL.track, fullLength, seqServer, 0, fullLength, GVB_GLOBAL.maxSize, mean, stdev, processedData.max, GVB_GLOBAL.ws, processedData.bins);

    if(processedData.hasOwnProperty("search") && processedData.search.points.hasOwnProperty(chromosome))
        {
        console.log("There's a search!")
        drawSearch(processedData.search);
        }
}




// SEARCH POINTS
////////////////////////////////
function searchPattern()
{
    if(globalDL1.drawn)
    {
        var pattern     = $('#patternSearch').val();
        var d           = $('#dSearch').val();

        var geo="none"; //can be gene, 5' UTR, 3' UTR, ncRNA gene, exons or intergene regions (neither gene neither intergene include pseudogenes or ncRNA genes)
        if(GVB_GLOBAL.geo)
            {
            var e = document.getElementById("geo_type");
            geo = e.options[e.selectedIndex].value;
            }

        if (DEBUG_GBV) console.log("\n----- SEARCH -----");
//        Server.search(drawSearch, pattern,d, geo, GVB_GLOBAL.intersect);
        Server.search(drawSearch, pattern,d, geo, "soft", GVB_GLOBAL.softMutations);
    }
}


// DRAW POINTS ON DATALINE 1
////////////////////////////////
function drawSearch(result)
{
    /**
     * @typedef {Object} result
     * @property points
     * @property sizePattern
     */

    // DATALINE 1: DRAW POINTS
    //----------------------------------
    dataLine_1_drawPoints(result.points, result.sizePattern);


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

    globalTime=new Date();

    if(typeof(sizePattern)=="object")//TODO: allow different sizes for each point
        sizePattern=30;
    Server.allAnnotationsGenes(getEnrichment, allPoints, "[\"gene\"]", sizePattern*ws, "left", "False", chromosomes, ws, "soft");
}


// GET ENRICHMENT
////////////////////////////////
function getEnrichment(gis, annotations)
{
    console.log("REAL TIME IN ANNOTATIONS: "+(new Date()-globalTime));
    setAnnotations(gis, annotations);

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
    if(DEBUG_GBV && clear)
    {
        // Reset input files
        $("#files").val('');

        console.log(new Array(15).join("\n"));
        console.log("Reset all (with File APIs)...");
    }
}
