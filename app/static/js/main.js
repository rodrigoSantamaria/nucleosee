/*
  ┌────────────────────────────────────────────────────────────┐
  │ main.js                                                    │
  ├────────────────────────────────────────────────────────────┤
  │ Description: Control logic for the browser                 │
  |  GPL v3 by Rodrigo Santamaría (University of Salamanca)    |
  └────────────────────────────────────────────────────────────┘
*/

var DEBUG_GBV = true;

var globalTime;

var GVB_GLOBAL =
{
    filenames: null,         // filenames: name of files
    files: null,
    chromosomes: null,      // chromosomes: array with the names of all chromosomes
    track: null,            // track: name of selected chromosome
    ws: 30,                 // window size: discrete to real ratio
    nb: 5,                  // num bins
    maxSize: 5000,        // maximum number of normalized data to store
    intersect: "soft",   // If true Hard instersects are considered (whole inclusion of pattern in genomic annotation for enrichment, etc.)
    geo: false,            // True if some genomic information (genes, UTRs, etc.) is used to filter out pattern matches
    softMutations: true,   // if true soft mutations (only 1-distance switch) is allowed to the patterns. E.g. "b" may change to "c" or "a" but not to "d" or "e".
    grid : false,           // if true, a grid to show percentile regions is shown on both lane 1 and 2
    forceReload : false,     //if true, available preprocessed data are discarded and the full preprocessing is remade based on actual parameters
    interpolation : "none",  //Either 'none', 'step' or 'rolling' for no interpolation, step-like or a rolling average 30-length window
    sendFile : false,    //if true, wig file should be uploaded to server for preprocessing
    exportFileURL: "",

    //Thesea are some variables that we should keep on track changes
    goURL: "",
    fastaURL: "",
    genesURL: "",
    posURL: "",
    gis: undefined,
    annotations: undefined
};


// CHECK FILE (AND UPLOAD FILE)
////////////////////////////////////////////

function checkFile(files)
{
    /**
     * @typedef {Object} file
     * @property name
     */

    destroyAll(false, false);

    // We ensure that there is a selected file
    if(typeof(files) != 'undefined')
    {
        GVB_GLOBAL.filenames=[];
        GVB_GLOBAL.files=[];

        for(var i=0;i<files.length;i++) {
            if (DEBUG_GBV) console.log("\n----- UPLOAD FILE -----");
            if (DEBUG_GBV) console.log("Name of file: " + files[i].name);


            GVB_GLOBAL.filenames.push(files[i].name);
            GVB_GLOBAL.files.push(files[i]);

            Server.testFile(files[i]);
            }
    }
}

/**
 * Asks the server for available preprocessed data
 */
function populateDataList(data)
    {
    var elSel = document.getElementById('selectionList');

    var cont=0
    var keys=Object.keys(data).sort(function (a, b) {
        return a.toLowerCase().localeCompare(b.toLowerCase());
        });

    for(var i in keys)
        {
        var d=data[keys[i]]
        var option = new Option(keys[i], keys[i]);
        option.title = "orig. file: "+d.file+"\norganism: "+d.organism+"\nwindow: "+ d.window+"\n#intervals:"+ d.intervals;
        elSel.options[elSel.options.length] = option;
        cont+=1;
        }

    elSel.size=Math.min(cont, 20);
    }

function populateOrganismList(data)
{
    var elSel = document.getElementById('speciesList');

    for(var i in data)
        elSel.options[elSel.options.length]=new Option(data[i],data[i]);
}

function selectData()
    {
    //0) Remove previous elements
    destroyAll(false, false);

    if($("#selectionList option:selected").length>=2)
        {
        $("#paramSearchDataset2")[0].disabled=false;
        $("#paramSearchDataset2a")[0].disabled=false;
        $("#paramSearchJoin")[0].disabled=false;
        $("#patternSearch2")[0].disabled=false;
        $("#paramAgnostic")[0].disabled=false;
        }

    //1) Get and print data
    for(var i=0; i<$("#selectionList option:selected").length; i++)
        {
        showImageLoading("imgLoadingFile", true);
        $('#loadText')[0].innerHTML="loading data...";

        //populate search option lists
        var elSel = document.getElementById('paramSearchDataset');
        var elOptNew = document.createElement('option');
        elOptNew.text = $("#selectionList option:selected")[i].value;
        elOptNew.value = $("#selectionList option:selected")[i].value;
        elSel.add(elOptNew, null);
        if(i==0)
            elSel.value=elOptNew.value;

        var elSel = document.getElementById('paramSearchDataset2');
        var elOptNew = document.createElement('option');
        elOptNew.text = $("#selectionList option:selected")[i].value;
        elOptNew.value = $("#selectionList option:selected")[i].value;
        elSel.add(elOptNew, null);
        if(i==1)
            elSel.value=elOptNew.value;

        //populate search option lists
        var elSel = document.getElementById('paramSearchDataset1a');
        var elOptNew = document.createElement('option');
        elOptNew.text = $("#selectionList option:selected")[i].value;
        elOptNew.value = $("#selectionList option:selected")[i].value;
        elSel.add(elOptNew, null);
        if(i==0)
            elSel.value=elOptNew.value;

        var elSel = document.getElementById('paramSearchDataset2a');
        var elOptNew = document.createElement('option');
        elOptNew.text = $("#selectionList option:selected")[i].value;
        elOptNew.value = $("#selectionList option:selected")[i].value;
        elSel.add(elOptNew, null);
        if(i==1)
            elSel.value=elOptNew.value;
            //elSel.selected(1);
        console.log("selection done");
        var clear="true";
        if(i>0)
            clear="false"
        Server.selectData(drawingFirstDataLine, $("#selectionList option:selected")[i].value, i, $("#selectionList option:selected").length, clear);
        //Server.getTrack(drawingFirstDataLine, GVB_GLOBAL.track,  $("#selectionList option:selected")[i].value, i, $("#selectionList option:selected").length);
        }

    }

function sendFile()
    {
    if (DEBUG_GBV) console.log("\n----- SENDING FILE -----");
    if(GVB_GLOBAL.sendFile)
        {
            $('#loadText')[0].innerHTML="sending file...";
            Server.sendFile(preprocessing, GVB_GLOBAL.files, 0)
        }
    }

// PREPROCESSING: NORMALIZE, STATISTICAL DESCRIPTORS, DISCRETIZE...
/////////////////////////////////////////////////////////////////////////
function preprocessing(chromosome)
{
    GVB_GLOBAL.sendFile=false;
    if (DEBUG_GBV) console.log("\n----- PREPROCESSING -----");
    if(typeof(chromosome) === 'undefined')
    {
        GVB_GLOBAL.track = "None";
        if (DEBUG_GBV) console.log("chromosome: " + GVB_GLOBAL.track + " (first chromosome found)");

        $('#loadText')[0].innerHTML="preprocessing data...";


        //populate search option lists
        var elSel = document.getElementById('paramSearchDataset');
        var elOptNew = document.createElement('option');

        elOptNew.text = $('#paramDescription').val();
        elOptNew.value = $('#paramDescription').val();
        elSel.add(elOptNew, null);

        var elSel = document.getElementById('paramSearchDataset2');
        var elOptNew = document.createElement('option');
        elOptNew.text = $('#paramDescription').val();
        elOptNew.value = $('#paramDescription').val();
        elSel.add(elOptNew, null);

        //populate search option lists
        var elSel = document.getElementById('paramSearchDataset1a');
        var elOptNew = document.createElement('option');

        elOptNew.text = $('#paramDescription').val();
        elOptNew.value = $('#paramDescription').val();
        elSel.add(elOptNew, null);

        var elSel = document.getElementById('paramSearchDataset2a');
        var elOptNew = document.createElement('option');
        elOptNew.text = $('#paramDescription').val();
        elOptNew.value = $('#paramDescription').val();
        elSel.add(elOptNew, null);

        Server.preprocess(filePreprocessed, GVB_GLOBAL.filenames, GVB_GLOBAL.track, $("#paramWS").val(), $("#paramNB").val(),
            GVB_GLOBAL.maxSize, $("#speciesList")[0][$("#speciesList")[0].selectedIndex].value,
            $("#interpolationList")[0][$("#interpolationList")[0].selectedIndex].value, $("#paramSD").val(),
            $("#paramDescription").val());
    }
    else    //track selection
    {
        GVB_GLOBAL.track = chromosome;
        if (DEBUG_GBV) console.log("chromosome: " + GVB_GLOBAL.track);

        //0) Remove previous elements
        destroyAll(false, true);

        //1) Get and print data
        for(var i=0; i<$("#selectionList option:selected").length; i++)
        {
            showImageLoading("imgLoadingFile", true);
            $('#loadText')[0].innerHTML="loading data...";
            Server.getTrack(drawingFirstDataLine, GVB_GLOBAL.track,  $("#selectionList option:selected")[i].value, i, $("#selectionList option:selected").length);
        }
    }
}

function filePreprocessed(message)
    {
    GVB_GLOBAL.sendFile=true;//for next files
    $('#loadText')[0].innerHTML=message;
    }


// DRAWING: DATALINE 1
////////////////////////////////
/**
 *
 * @param processedData - data to draw
 * @param chromosome - track (usually chromosome) to be drawn in foreground
 * @param index - in de case that several processedData are to be loaded, which one is this one
 */
function drawingFirstDataLine(processedData, chromosome, index, total) {
    /**
     * @typedef {Object} processedData
     * @property seq
     * @property mean
     * @property stdev
     * @property dseq
     * @property fullLength
     * @property chromosomes
     */
    if(index>=total-1) {
        showImageLoading("imgLoadingFile", false);
        $('#loadText')[0].innerHTML = "Data loaded";
        }


    GVB_GLOBAL.forceReload = false;

    // We create icons chromosomes and bind the click event
    GVB_GLOBAL.chromosomes = processedData.chromosomes;
    if (chromosome == "None")
        {
        createIconsChromosomes(GVB_GLOBAL.chromosomes);
        GVB_GLOBAL.track = processedData.chromosomes[0];
        }
    else
        GVB_GLOBAL.track = chromosome;


    var seqServer = processedData.seq;    // just a sampling of about 400K of the original full length sequence
    var fullLength = processedData.fullLength;

    // DATALINE 1 (preprocessed data)
    //----------------------------------
    if (DEBUG_GBV) console.log("\n----- DRAWING: DATALINE 1 (" + GVB_GLOBAL.track + ")-----");
    if (DEBUG_GBV) console.log("Length of seqServer:" + seqServer.length + " (full length seq=" + fullLength + ")");


    if (typeof(fullLength) == "object")
        dataLine_1(processedData, 0, fullLength[GVB_GLOBAL.track], total, GVB_GLOBAL.annotations, GVB_GLOBAL.gis);
    else
        dataLine_1(processedData, 0, fullLength, total, GVB_GLOBAL.annotations, GVB_GLOBAL.gis);

    //In case of switching tracks
    if (processedData.hasOwnProperty("search") && processedData.search[globalSeqs[index].dataName].points.hasOwnProperty(chromosome))
        drawSearch(processedData.search);
    if (processedData.hasOwnProperty("ego"))
        drawEnrichment(processedData.ego);
}




// SEARCH POINTS
////////////////////////////////
/**
 * Search workflow is as follows:
 *
 *          pattern --> SEARCH --> positions --> GET_ALL_ANNOT_GENES --> gids --> ENRICHMENT
 * API                  search                   annotations                      enrichmentGO
 * SOURCE               wig                      gff                              go+goa
 *
 *
 *  There's an additional workflow on navigation on positions
 *              single_position --> ANNOTATIONS_GENES  --> gene details
 *  CLIENT                          draw.dataLine_2
 *  API                             annotations
 *  SOURCE                          gff
 *
 *
 * And also for sequences. There we can try to outsurce to eutils (example with S cerevisiae):
 *
 * 1) Get organism id (in id_list)
 * http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term=Saccharomyces%20cerevisiae&retmode=json  --> note this will provide for ALL strains of Sc
 *
 * 2) Get sequence id (parse query key and webenv)
 * http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=genome&db=nuccore&id=15&term=chromosome&cmd=neighbor_history  --> we can use just a taxon id
 * http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=genome&db=nuccore&id=4392&term=chromosome&cmd=neighbor_history --> it still generates the whole Sc family fasta! (1.5G)
 *
 * 3) Get the real stuff
 * http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key=2&WebEnv=NCID_1_44089431_130.14.22.215_9001_1469092699_175516958_0MetA0_S_MegaStore_F_1&rettype=fasta&retmode=txt
 *
 * id can be 0 if it comes from the quick search or 1 if it comes from the advanced search
 */
function searchPattern(id)
{
    var selection1=$("#paramSearchDataset option:selected")[0].value;
    var selection2=$("#paramSearchDataset2 option:selected")[0].value;
    var selectionJoin=$("#paramSearchJoin option:selected")[0].value;

    var globalDL1=getLevel(selection1, "DL1");

    if(globalDL1.drawn)
    {

        if(id==1)
            var pattern = $('#patternSearch').val();
        else    //simple search (no join), use d as in the complicate search
            {
            var pattern = $('#patternSearch0').val();
            selectionJoin="None"
            }
        var pattern2 = $('#patternSearch2').val();
        var d = $('#dSearch').val();

        var geo = "none"; //can be gene, 5' UTR, 3' UTR, ncRNA gene, exons or intergene regions (neither gene neither intergene include pseudogenes or ncRNA genes)
        if (document.getElementById('paramGeo').checked) {
            var e = document.getElementById("geo_type");
            geo = e.options[e.selectedIndex].value;
        }
        var intersect = "soft";
        if (document.getElementById('paramIntersect').checked)
            intersect = "hard"//TODO: hard with intergenic regions is not working.

        var agnostic="false";
        var section=300;
        var portion=150;
        if(document.getElementById('paramAgnostic').checked)
            {
            agnostic = "true";
            section=$('#agnosticSearchL').val();
            portion=$('#agnosticSearchV').val();
            selection1=$("#paramSearchDataset1a option:selected")[0].value;
            selection2=$("#paramSearchDataset2a option:selected")[0].value;
            }

        if (DEBUG_GBV) console.log("\n----- SEARCH -----");
        Server.search(searchResults, pattern, d, geo, intersect, GVB_GLOBAL.softMutations,
            selection1, selection2, selectionJoin, pattern2,
            agnostic, section, portion);
    }
}


// DRAW POINTS ON DATALINE 1
////////////////////////////////
function searchResults(result)
{
    /**
     * @typedef {Object} result
     * @property points
     * @property sizePattern
     */

    // DATALINE 1: DRAW POINTS
    //----------------------------------
    //drawSearch(result.points, result.sizePattern, dataName);
    drawSearch(result);

    exportPositions(result)
    //Server.exportFASTA(exportFasta, globalSeqs[0].dataName);

    // ENRICHMENT
    //----------------------------------
    if (DEBUG_GBV) console.log("\n----- ENRICHMENT -----");
    getAllAnnotations(result);
}

function exportFasta(text)
    {
    var blob = new Blob([text], {type: 'text/plain'});
    var URL = this.URL || this.webkitURL;
    var textURL;
    if (textURL !== null)   URL.revokeObjectURL(textURL);
    var textURL = URL.createObjectURL(blob);
    dl1[0].fastaURL=textURL;
    window.open(dl1[0].fastaURL)
    }

/**
 * Export positions in BED format: chromosome start end
 * @param points
 */
function exportPositions(points)
    {
    //var text="pos\tchromosome\tdataset\n"; //no header for bed format
    var text="";
    for(var dataName in points)
        {
        var size=points[dataName]["sizePattern"]
        size=size*globalSeqs[0].ws;
        for(var ch in points[dataName]["points"])
            {
            var pch = JSON.parse(points[dataName]["points"][ch]);
            for (var p in pch)
                {
                //text = text +  ch + "\t" + pch[p] + "\t" + dataName +"\n";
                text = text +  ch + "\t" + pch[p] + "\t"+ (pch[p]+size) +"\n";
                }
            }
        break;//we only pick for one (maybe buggy, fix this)
        }
    var data = new Blob([text], {type: 'text/plain'});
    var URL = this.URL || this.webkitURL;
    var textURL;
    if (textURL !== null)   URL.revokeObjectURL(textURL);

    var textURL = URL.createObjectURL(data);
    var globalDL1=dl1[0];

    //Clean up previous urls
    globalDL1.genesURL="http://vis.usal.es/seqview/client/html/error/404.html";
    globalDL1.goURL="http://vis.usal.es/seqview/client/html/error/404.html";
    globalDL1.fastaURL="http://vis.usal.es/seqview/client/html/error/404.html";
    globalDL1.posURL=textURL;
    }

function exportGenes(annotations)
{
    var text="id\tpos\tstart\tend\tsense\tchromosome\tname\n";
    for(var a in annotations)
    {
        var aa=annotations[a]
        text=text+a+"\t"+aa.pos+"\t"+aa.start+"\t"+aa.end+"\t"+aa.sense+"\t"+aa.chromosome+"\t"+aa.name+"\n";
    }
    var data = new Blob([text], {type: 'text/plain'});
    var URL = this.URL || this.webkitURL;
    var textURL;
    // If we are replacing a previously generated file we need to
    // manually revoke the object URL to avoid memory leaks.
    if (textURL !== null) {
        URL.revokeObjectURL(textURL);
    }

    var textURL = URL.createObjectURL(data);

    var globalDL1=dl1[0];//Over the first of all the level1 lines (usually just one)
    globalDL1.genesURL=textURL;
    }


function exportGO(annotations)
{
    var text="id\tname\t#gis\t#go\tp-value\tgis\n";
    for(var a in annotations)
    {
        var aa=annotations[a]
        text=text+a+"\t"+aa.go_name+"\t"+aa.ngis+"\t"+aa.ngo+"\t"+aa.pval+"\t"+aa.gis+"\n";
    }
    var data = new Blob([text], {type: 'text/plain'});
    var URL = this.URL || this.webkitURL;
    var textURL;
    // If we are replacing a previously generated file we need to
    // manually revoke the object URL to avoid memory leaks.
    if (textURL !== null) {
        URL.revokeObjectURL(textURL);
    }

    var textURL = URL.createObjectURL(data);
    var globalDL1=dl1[0];
    globalDL1.goURL=textURL;
}

function exportFASTA()
    {
        var text="";
        for(var a in annotations)
        {
            var aa=annotations[a]
            text=text+a+"\t"+aa.go_name+"\t"+aa.ngis+"\t"+aa.ngo+"\t"+aa.pval+"\t"+aa.gis+"\n";
        }
        var data = new Blob([text], {type: 'text/plain'});
        var URL = this.URL || this.webkitURL;
        var textURL;
        // If we are replacing a previously generated file we need to
        // manually revoke the object URL to avoid memory leaks.
        if (textURL !== null) {
            URL.revokeObjectURL(textURL);
        }

        var textURL = URL.createObjectURL(data);
        var globalDL1=dl1[0];
        globalDL1.goURL=textURL;
    }

// GET ALL ANNOTATIONS
////////////////////////////////
function getAllAnnotations(matches)
{
    var chromosomes = GVB_GLOBAL.chromosomes;
    var ws          = GVB_GLOBAL.ws;

    var extend=parseInt($("#paramExtend").val()); //this parameter is to allow extend the area to search for GIS x nucleotides on left and right from the pattern size searched

    globalTime=new Date();

    var allPoints={}
    var sizePattern=undefined;

    for(var i in matches)
        if(sizePattern==undefined)
            sizePattern=matches[i].sizePattern
        for(var j in matches[i].points) {
            var pos=matches[i].points[j].replace("[","").replace("]", "").split(",")
            if(pos[0]=="")
                pos=[]

            for(var k in pos)
                pos[k]=parseInt(pos[k])-extend
            //var pos=matches[i].points[j]-extend;//extend on the left
            if (allPoints[j] == undefined)
                //allPoints[j] = matches[i].points[j];
                allPoints[j] = pos;
            else
                //allPoints[j].push(matches[i].points[j]);
                allPoints[j].push(pos);
            }


    var winS;
    if(typeof(sizePattern)=="object")//TODO: allow different sizes for each point
        {
        winS={};
        for(var i in sizePattern)
            {
            winS[i]=JSON.parse(sizePattern[i]);
            for(var j in winS[i])
                {
                winS[i][j] *= ws;
                winS[i][j] +=extend; //extend on the left
                }
            winS[i]="["+winS[i].join(",")+"]";
            }
        }
    else
        winS=sizePattern*ws;
    Server.allAnnotationsGenes(getEnrichment, allPoints, "[\"any\"]", winS, "left", "False", chromosomes, ws, "soft", globalSeqs[0].dataName, globalSeqs[0].track);
}


// GET ENRICHMENT
////////////////////////////////
function getEnrichment(gis, annotations)
{
    console.log("REAL TIME IN ANNOTATIONS: "+(new Date()-globalTime));
    setAnnotations(gis, annotations);

    //Prepare file for export
    exportGenes(annotations);

    //Server.enrichmentGO(drawEnrichment, gis, "fdr", 0.01)

    var pval = $('#paramPval').val();

    var discard="[]"
    if(document.getElementById('paramDiscard').checked)
        discard="[\"IEA\"]";


    Server.enrichmentGO(drawEnrichment, "", "fdr", pval, discard)
}








function createIconsChromosomes(chromosomes)
{
    $('#imagesChromosomes').empty();
    for(var i=0; i<chromosomes.length; i++)
    {
        $('#imagesChromosomes').append('<img style="margin-top:15px;margin-right:5px"'+
            'id="'+chromosomes[i]+'" class="image-chromosome" data-chromosome="'+chromosomes[i]+'" '+
            'src="images/chromosome.png" height="24px" width="24px">');
            //TODO: ponerles un tamaño proporcional a su longitud puede molar
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

function destroyAll(clear, tracksOnly)
{
    // Empty all SVG images (of the array)
    if(dl1!=undefined)
        for(var i=0;i<dl1.length;i++)
            $("#"+dl1[i].cv.nameSVG).empty();

    if(dl2!=undefined)
        for(var i=0;i<dl2.length;i++)
            $("#"+dl2[i].cv.nameSVG).empty();
    var images = ["lineSeq3"];
    for(var i=0;i<images.length;i++)
    {
        var image = $("#"+images[i]);
        if(image.html() != "")
            image.empty();
    }



    // Empty all icons of chromosomes
    if(tracksOnly==false) {

        $("#imagesChromosomes").empty();

        // Clear console and files configuration
        if (DEBUG_GBV && clear) {
            // Reset input files
            $("#files").val('');

            console.log(new Array(15).join("\n"));
            console.log("Reset all (with File APIs)...");
            }

        //Clear data pattern search selection lists
        document.getElementById('paramSearchDataset').options.length=0;
        document.getElementById('paramSearchDataset2').options.length=0;
        document.getElementById('paramSearchDataset1a').options.length=0;
        document.getElementById('paramSearchDataset2a').options.length=0;

        GVB_GLOBAL.annotations=undefined
        GVB_GLOBAL.gis=undefined
    }
    else    //save annotations for further usage
        {
        GVB_GLOBAL.annotations=dl1[0].annotations
        GVB_GLOBAL.gis=dl1[0].gis
        GVB_GLOBAL.goURL=dl1[0].goURL
        GVB_GLOBAL.posURL=dl1[0].posURL
        GVB_GLOBAL.genesURL=dl1[0].genesURL
        GVB_GLOBAL.fastaURL=dl1[0].fastaURL
        }

    globalSeqs=[]
    dl1=[]
    dl2=[]
}

