/*
 ┌────────────────────────────────────────────────────────────┐
 │ draw.js                                                    │
 ├────────────────────────────────────────────────────────────┤
 │ Description: View logic for the browser (frontend)         │
 |  GPL v3 by Rodrigo Santamaría (University of Salamanca)    |
 └────────────────────────────────────────────────────────────┘
 */



var marginDL =
{
    top:10,
    bottom:60,
    left:50,
    right:20
};

var dimDL =
{
    graphWidth : screen.width,
    graphHeight : 150,

    width : screen.width-marginDL.right-marginDL.left,
    height : 150-marginDL.top-marginDL.bottom
};

var dimDLAnnotation =
{
    x0:5,
    y0:7,
    lineHeight : 12
};

var DLcolors=["steelblue", "#FF4400", "#968096", "#327A2A", "#70AA00", "#FAB526"];
var DLcolors=["steelblue", "#968096", "#327A2A", "#70AA00", "#FAB526"];
//var DLcolors=["steelblue", "#327A2A", "#968096", "#70AA00", "#FAB526"];
var currentZoom=1;

//________________________ SEQUENCE OBJECTS DEFINITION _____________________________

//This object compiles the summay information for a given wig sequence
//globalDL1,2 and 3 contain the information of different detail levels
function globalSeqObject() {
    this.scaleSeqServ = null;//Scale between sequence and server
    this.seqServ = null;//Sequence for current track as provided by the server (sampled to screen width)
    this.mean = null;//Mean value for current track
    this.stdev = null;//Standard deviation for current track
    this.max = null;//Maximum value for current track
    this.min = null;//Minimum value for current track

    //The following parameters should be the same for all the sequence to be loaded
    this.fullLength = null;//number of nucleotides for current track
    this.bins = null;//number of bins for BWT discretization
    this.ws = null;//window size for BWT discretization
    this.tracks = null;//all track names
    this.track = null;//current track

    this.dataName = null; //Dataset name to which this sequence belongs

    this.numReplicates = null;
    this.correlations = null;
    this.id=null;
    }
var globalSeqs=[]
//In drawSearch, drawEnrichment, and dataLine2, any sequence will do (ref to globalSeqs[0]) as they only use ws, bins, track(s)


function marginObject()
    {
    this.top=marginDL.top;
    this.bottom=marginDL.bottom; //bottom should be smaller if several DL1 tracks are drawn
    this.left=marginDL.left;
    this.right=marginDL.right;
    }

function dimObject()
    {
        this.graphWidth = dimDL.graphWidth;
        this.graphHeight = dimDL.graphHeight;
        this.width = dimDL.width;
        this.height = dimDL.height;
    }
function dimAnnotationObject()
    {
        this.x0 = dimDLAnnotation.x0;
        this.y0 = dimDLAnnotation.y0;
        this.lineHeight =dimDLAnnotation.lineHeight;
    }
function customViewObject(name, className, bracketHeight)
    {
        this.nameSVG = name;    //SVG name for the HTML
        this.classSVG = className;       //class name for the CSS

        this.svg = null;             //svg object with the image
        this.scaleSeqScreen = null;  //scale ratio between compressed sequence and screen width
        this.scaleServScreen = null; //scale ratio between full (server) sequence and screen width
        this.data = null;            //
        this.xScale = null;          //scaling function for X
        this.yScale = null;          //scaling function for Y
        this.color="steelblue";

        this.margin = new marginObject();          //pixels on each side
        this.dim = new dimObject();              //graph dimensions
        this.dimAnnotation = new dimAnnotationObject();    //annotations dimensions (for enriched GO terms)

        this.bracketHeight=bracketHeight; //actually only used in DL2
    }

//Saves information about a set of sequences
function sequencesObject()
{
    this.seqs=null; //actual nucleotide sequences
    this.alignment=null;    //alignment proposed (currently not used as it's computationally costly)
    this.aconsensus=null;   //consensus sequence
    this.aprofile=null; //consensus profile
    this.consensus=null;
    this.profile=null;

    this.motifs=null;   //motifs matching the consensus alignment
    this.locations=null;    //locations of the motifs on the sequence
    this.motifConsensus=null;   //consensus of all the motifs
    this.motifProfile=null; //
}


// Chromosome level. Level 1 saves information about the chromosome level, mostly visual information in cv
// (margins, scales, init point) and metainfo for searches (match points, gis, annotations)
function globalDL1Object(id)
{
    // Core variables (used in the core)
    this.cv= new customViewObject(id, "dl1", null);
    this.drawn = false;

    //this.seqPoints = null;   //In the case of a search, the array of matching points (UNUSED)
    this.pointTip = null;

    this.gis = null;         //In the case of a search, the genes that include matching positions
    this.annotations = {};   //In the case of a search, the enriched annotation for the genes that include matching positions
    this.genesURL= "http://vis.usal.es/seqview/client/html/error/404.html"; //URLs for raw text to matching positions, genes and GO terms.
    this.posURL= "http://vis.usal.es/seqview/client/html/error/404.html";
    this.goURL= "http://vis.usal.es/seqview/client/html/error/404.html";
    this.fastaURL= "http://vis.usal.es/seqview/client/html/error/404.html";
}


// Nucleosome level (level 2) saves information at nucleosome level, usualy 1:1 scale respect to real data
function globalDL2Object(id)
{

    // Core variables (used in the core)
    this.cv= new customViewObject(id, "dl2", 40);
    this.drawn=false;

    this.startSeq=null;//starting and ending positions respect to real data (defines the interval)
    this.endSeq=null;
};


// Nucleotide level (level 3) saves information about the actual sequences as retrieved from FASTA files. Should
//consider enough pixel width to make letters readable
function globalDL3Object()
{
    // Core variables (used to draw)
    this.cv= new customViewObject("lineSeq3", "dl3", 120);
    this.sequences=new sequencesObject();
};

var dl1=[];//array of dl1 elements
var dl2=[];//array of dl2 elements
var globalDL3=new globalDL3Object();



//-------------------------------------------------------------
//          DATALINE 1
//-------------------------------------------------------------
/**
 *
 * @param processedData data to draw
 * @param startSeq  starting point
 * @param endSeq    ending point
 * @param numLines  total number of tracks for this level (total number of different data, usually 1 or 2)
 */
function dataLine_1(processedData, startSeq, endSeq, numLines, annot, gis)
{
    var startTime = new Date();
    var globalSeq = new globalSeqObject();

    var fullLength=processedData.fullLength[GVB_GLOBAL.track];

    // Save information of the genomic sequence
    globalSeq.tracks        = processedData.chromosomes; // all chromosomes
    globalSeq.track         = GVB_GLOBAL.track; // chromosome or wig track it refers to
    globalSeq.seqServ       = processedData.seq;
    globalSeq.mean          = processedData.mean;
    globalSeq.stdev         = processedData.stdev;
    globalSeq.max           = processedData.max;
    globalSeq.min           = processedData.min;
    globalSeq.bins          = processedData.bins;
    globalSeq.dataName      = decodeURIComponent(processedData.dataName);
   // globalSeq.ws            = GVB_GLOBAL.ws;
    globalSeq.ws            = processedData.ws;
    globalSeq.scaleSeqServ  = 1;
    globalSeq.id            = globalSeqs.length;
    if(Math.floor(fullLength/GVB_GLOBAL.maxSize) >= 1)
        globalSeq.scaleSeqServ = Math.floor(fullLength/GVB_GLOBAL.maxSize);

    globalSeqs.push(globalSeq);


    var globalDL1=new globalDL1Object("lineSeq_"+dl1.length);
    dl1.push(globalDL1);

    var globalDL2=new globalDL2Object("lineSeq2_"+dl2.length);
    dl2.push(globalDL2);
    if(dl1.length==1 && annot!=undefined)
        setAnnotations(gis, annot);
    if(dl1.length==1 && GVB_GLOBAL.posURL!="")
        {
        dl1[0].posURL=GVB_GLOBAL.posURL;
        dl1[0].genesURL=GVB_GLOBAL.genesURL;
        dl1[0].goURL=GVB_GLOBAL.goURL;
        dl1[0].fastaURL=GVB_GLOBAL.fastaURL;
        }

    if(dl1.length>1)
        {
        dl1[dl1.length-1].cv.color=DLcolors[dl1.length-1];
        dl2[dl2.length-1].cv.color=DLcolors[dl2.length-1];
        globalDL1=dl1[dl1.length-1];
        globalDL2=dl2[dl2.length-1];
        }

    // We use the core function
    if(dl1.length==1)
        dataLine_core(globalDL1, globalSeq, globalSeq.seqServ, globalSeq.scaleSeqServ, startSeq, endSeq, startSeq);
    else
        dataLine_lineOnly(globalDL1, globalSeq, globalSeq.seqServ, globalSeq.scaleSeqServ, startSeq, endSeq, startSeq);


    // We confirm that we have finished
    globalDL1.drawn = true;

    if(DEBUG_GBV) console.log("Time spent dataLine1: "+ (new Date()-startTime)+"ms");
}



function getLevel(dataName, level, realDL1)
    {
    for(var i in globalSeqs)
        if(globalSeqs[i].dataName==dataName) {
            if (level == "DL1")
                if(realDL1==undefined)
                    return dl1[0]
                else
                    return dl1[i]
            if (level == "DL2")
                return dl2[i]
            if (level == "Seq")
                return globalSeqs[i]
            }
    return null;
    }



function drawSearch(allPoints0)
{
    var matches = [];

    //remove all previous points
    for(var i in dl1)
        if(dl1[i].cv.svg!=undefined)
             dl1[i].cv.svg.selectAll(".point").remove();

    // Remove all occurrences label
    dl1[0].cv.svg.selectAll(".search-label").remove();
    dl1[0].cv.svg.selectAll(".export-label").remove();

    //For different possible sizes in different datasets
    for(var i in Object.keys(dl2))
        dl2[i].sizePattern=0;
    for(var dataName in allPoints0)
        {
        var globalDL2=getLevel(dataName, "DL2")
        globalDL2.sizePattern=allPoints0[dataName].sizePattern
        }

    var dataName=Object.keys(allPoints0)[0];
    //for(var dataName in allPoints0)
        {//FOR START
        var allPoints=allPoints0[dataName].points;
        var sizePattern=allPoints0[dataName].sizePattern;
        var globalSeq = getLevel(dataName, "Seq");
        var globalDL1 = getLevel(dataName, "DL1");


        var track = globalSeq.track;
        var tracks = globalSeq.tracks;
        var points = JSON.parse(allPoints[track]);
        var numNucleotidesDraw = globalDL1.cv.dim.width; // because the scale is 1:1
        if (typeof(sizePattern) == "object")
            sizePattern = JSON.parse(sizePattern[track]);
        if (DEBUG_GBV) console.log("\ndataLine_1_drawPoints(): " + points.length + " points");

        for (i in tracks) {
            var pp = JSON.parse(allPoints[tracks[i]]);
            if(matches[i]==undefined)   matches[i]=pp.length;
            else                        matches[i]+=pp.length;
        }

        // We calculate the points found in the sequence
        var seqPoints = [];
        for (var i = 0; i < points.length; i++)
            seqPoints.push(points[i]);

        // Now, with seqPoints, we calculate the points found in the dataLine
        var dataPoints = [];
        for (i = 0; i < seqPoints.length; i++) {
            var dataPoint = Math.round(seqPoints[i] / globalDL1.cv.scaleSeqScreen);
            if (dataPoint >= globalDL1.cv.data.length)
                console.log("Out of bounds!: " + dataPoint + ", " + seqPoints[i]);
            else {
                if (typeof(sizePattern) == "number")//in the case of BWT/interval searches
                    dataPoints.push({
                        real_pos: seqPoints[i],
                        pos: (globalDL1.cv.data)[dataPoint].pos,
                        value: (globalDL1.cv.data)[dataPoint].value,
                        size: sizePattern,
                        selected: false
                    });
                else                //in the case of nominal gene/go searches
                    dataPoints.push({
                        real_pos: seqPoints[i],
                        pos: (globalDL1.cv.data)[dataPoint].pos,
                        value: (globalDL1.cv.data)[dataPoint].value,
                        size: sizePattern[i],
                        selected: false
                    });
            }
        }

        // Mouseover tip and drawing the corresponding line
        globalDL1.pointTip = d3.tip()
            .attr('class', globalDL1.cv.classSVG + ' point-tip')
            .offset([30, 0]) // [top, left] to center the tip
            .html(function (d) {
                var point = d.real_pos;


                // We round to the nearest hundred from sizePattern. E.g. 270 --> 300,  540 --> 500
                var focusLine = Math.round(d.size * globalSeq.ws / 100) * 100;
                if (focusLine > dl2[0].cv.dim.width)
                    numNucleotidesDraw = focusLine * 1.1;//make it a 10% larger than the thing to draw
                else
                    numNucleotidesDraw = dl2[0].cv.dim.width;  //1:1 scale in the rest of cases
                var startSeq = point - (numNucleotidesDraw / 2) + (focusLine / 2);
                var endSeq = point + (numNucleotidesDraw / 2) + (focusLine / 2);

                //Remove nucleotide track
                removeNucleotides();

                zoom(currentZoom, startSeq, endSeq);

                return "<strong>" + d3.format(",")(point) + ":</strong> " + d3.format(".2f")(d.value); // e.g. "307,770 : 0.79"
            });

        // Calls tip
        globalDL1.cv.svg.call(globalDL1.pointTip);

        //Draw the points
        drawPoints(dataPoints, globalDL1);
    }//FOR EACH LANE


    // Draw occurrences label (only on first lane) --------------------
    var startX = 210;

    // In the case of Firefox, we require a bit more offset
    var ua= navigator.userAgent, tem,
    M= ua.match(/(opera|chrome|safari|firefox|msie|trident(?=\/))\/?\s*(\d+)/i) || [];
    if(/trident/i.test(M[1])){
        tem=  /\brv[ :]+(\d+)/g.exec(ua) || [];
        return 'IE '+(tem[1] || '');
    }
    if(M[1]=== 'Firefox'){
        startX=startX+70;
    }

    dl1[0].cv.svg.append("g")
        .attr("class", dl1[0].cv.classSVG + " export-label")
        .append("text")
        .attr('x', startX)
        .attr('y', -2)
        .on("click", function () {
            window.open(dl1[0].posURL);
        })
        .text("matches>");

    dl1[0].cv.svg.selectAll(".searchLabels")
        .data(matches)
        .enter()
        .append("text")
        .attr("class", dl1[0].cv.classSVG + " search-label")
        .attr('x', function (d, i) {
            return startX + 55 + i * 29;
        })
        .attr('y', -2)
        .style("font-weight", function (d, i) {
            if (tracks[i] == track) return "bold";
            return "";
        })
        .text(function (d) {
            return d + ""
        });


    dl1[0].cv.svg.append("g")
        .append("text")
        .attr("class", dl1[0].cv.classSVG + " export-label")
        .attr('x', startX + 55 + 29 * (matches.length))
        .attr('y', -2)
        .on("click", function () {
            window.open(dl1[0].genesURL);
        })
        .text("genes>");

    dl1[0].cv.svg.append("g")
        .append("text")
        .attr("class", dl1[0].cv.classSVG + " export-label")
        .attr('x', startX + 55 + 29 * (matches.length) + 40)
        .attr('y', -2)
        .on("click", function () {
            window.open(dl1[0].goURL);
        })
        .text("go>");

    dl1[0].cv.svg.append("g")
        .append("text")
        .attr("class", dl1[0].cv.classSVG + " export-label")
        .attr('x', startX + 55 + 29 * (matches.length) + 40 + 20)
        .attr('y', -2)
        .on("click", function () {
           if(dl1[0].fastaURL=="http://vis.usal.es/seqview/client/html/error/404.html")
                Server.exportFASTA(exportFasta, globalSeqs[0].dataName);
            else
                window.open(dl1[0].fastaURL)
        })
        .text("fasta>");

    //GET SEQUENCES, MOTIFS, (ALIGNMENT), CONSENSUS
    //NOTE: alignment takes more than 1s if there's >50 sequences! (using the fastest method: kalign) -> by now we disable it
    //NOTE: by now searches in seq are only possible with equal lengths
    //if (typeof(sizePattern) == "number" && points.length > 0)  //only in the case that all the searches have the same size (not for gene searches, etc. as align is costly
    //    Server.nucProfile(saveSequences, globalSeq.track, "[" + points + "]", sizePattern * globalSeq.ws, document.getElementById("paramKmotif").value, globalSeq.dataName);
    //else
    //    saveSequences(null);
}


function drawPoints(dataPoints, globalDL1)
    {
    //--- testing
    globalDL1=dl1[0];//No different places now
    //---
    globalDL1.cv.svg.selectAll("."+globalDL1.cv.classSVG+".point").remove();

    // Create all new points
    globalDL1.cv.svg.selectAll(".data")
        .data(dataPoints)
        .enter()
        .append("circle")
        .attr('class', globalDL1.cv.classSVG+' point')
        .attr('cx', function(d) { return globalDL1.cv.xScale(d.pos); })
        .attr('cy', function(d) { return globalDL1.cv.yScale(d.value); })
        .attr('r', function(d) {  return d.selected?3:2;      })
        .attr('stroke-width', function(d) {return d.selected?2:1;})
        .on('mouseover', globalDL1.pointTip.show)
        .on('mouseout', globalDL1.pointTip.hide);


    // Keep colored the selected point
    $("."+globalDL1.cv.classSVG+".point").bind( "mouseover", function()
    {
        $("."+globalDL1.cv.classSVG+".point").attr('class', globalDL1.cv.classSVG+' point');
        $(this).attr('class', globalDL1.cv.classSVG+' point pressed');
    });
    }

/**
 * Draws the names of enriched GO terms for the genes that contains the searched positions.
 * The terms are drawn with a size (inverselly) proportional to their p-value, and they can be hovered and clicked
 * @param enrichment
 */
function drawEnrichment(enrichment)
{
    // Remove all goterm (previous)
    //var globalDL1=dl1[dl1.length-1];//not in the latest of the DL1 lines
    var globalDL1=dl1[0];//not in the latest of the DL1 lines
    globalDL1.cv.svg.selectAll("."+globalDL1.cv.classSVG+".goterm").remove();
    globalDL1.cv.svg.selectAll("."+globalDL1.cv.classSVG+".goterm-tip").remove();
    d3.selectAll("."+globalDL1.cv.classSVG+".goterm-tip").remove();//some goterm-tips are not removed, not sure why

    //Return if no enrichment
    if(Object.getOwnPropertyNames(enrichment).length==0)
        return;

    //Prepare file for export
    exportGO(enrichment);


    // Adjust the scale
    var goSize = d3.scale.log().base(10)
        .domain([1,10e-6])  // p-values
        .clamp(true)
        .rangeRound([8, 16]); // min and max letter size
    //.rangeRounls
    // d([6, 14]); // min and max letter size

    // We create an array with goterms
    var goterms = [];
    //var keys=enrichment.keys().sort();
    for(var k in enrichment)
    //for(var k in keys)
    {
        if(enrichment.hasOwnProperty(k))
        {
            if (enrichment[k].go_name != undefined)
            {
                var e = enrichment[k];
                e["go_id"] = k;
                goterms.push(e);
            }
            else
                console.log("drawEnrichment(): Error! => GO term " + k + " not found (possibly outdated OBO file?");
        }

    }
    if(DEBUG_GBV) console.log("drawEnrichment(): "+goterms.length+" enriched terms");


    // Mouseover tip and show the information of goterm
    var tip = d3.tip()
        .attr('class', globalDL1.cv.classSVG+' goterm-tip')
        .offset([60, 0])   // [top, left] to center the tip
        .html(function(d,i)
        {
            var text = "p-value: "+ d3.format(".2e")(d.pval)+"<br>"+
                d.ngis+"/"+ d.ngo +" genes"+"<br>";
            return text;
        });

    // Calls tip
    globalDL1.cv.svg.call(tip);


    // Dimensions of annotations
    var x0 = globalDL1.cv.dimAnnotation.x0;
    var y0 = globalDL1.cv.dimAnnotation.y0;
    var lineHeight = globalDL1.cv.dimAnnotation.lineHeight;

    var dx = x0;
    var dy = [];
    var yline = 0;

    // Draw all goterm
    globalDL1.cv.svg.selectAll("."+globalDL1.cv.classSVG+".goterm")
        .data(goterms)
        .enter()
        .append("text")
        .attr('class', globalDL1.cv.classSVG+' goterm')
        .attr('x', function(d)
        {
            // We calculate the width of the text that we will write
            var widthText = getTextWidth(" "+d.go_name+" ·", goSize(d.pval)+"px sans-serif");

            // We found 'dx' and 'dy'
            dx += widthText;
            if(dx >= globalDL1.cv.dim.width)
            {
                yline++;
                dx = x0 + widthText;
            }
            dy[d.go_name] = y0 + globalDL1.cv.dim.height + 25 + lineHeight*yline;

            return dx-widthText;
        })
        .attr('y', function(d)
        {
            return dy[d.go_name];
        })
        .attr('font-size', function(d)
        {
            return goSize(d.pval)+"px";
        })
        .text(function(d,i)
        {
            // Name of goterm
            if(i != (goterms.length-1))
                return " "+d.go_name+" ·";
            else
                return " "+d.go_name;
        })
        .on('click', function(d)
        {

            //Underline the selected term
            //0) first remove previous selection
            //for(var i in dl1)
                {
              //  var gdl = dl1[i];
                var gdl=dl1[0];
                if (d3.select(this).classed("selected") == false)
                    gdl.cv.svg.selectAll("." + globalDL1.cv.classSVG + ".goterm").classed("selected", false);
                }
            d3.select(this).classed("selected")?d3.select(this).classed("selected", false):d3.select(this).classed("selected", true);

            //1) Determine points selected
            var gis=[];
            for(var i in d["gis"])
                {
                //var annot=globalDL1.annotations[d["gis"][i]];
                var annot=dl1[0].annotations[d["gis"][i]];//annotations are only in the first DL1 track
                var globalSeq=globalSeqs[0];
                if (annot["chromosome"] == globalSeq.track)
                    gis.push(annot["pos"]);
                }

            //2) Hightlight it
            //for(var k in dl1)
                {
                //var gdl=dl1[k];
                var gdl=dl1[0];
                if(gdl.cv.svg!=undefined) {
                    var dataPoints = gdl.cv.svg.selectAll("." + gdl.cv.classSVG + ".point").data();
                    if (d3.select(this).classed("selected") == true) {
                        for (var i in dataPoints) {
                            dataPoints[i].selected = false;
                            for (var j in gis) {
                                if (gis[j] > dataPoints[i]["real_pos"] - 50 && gis[j] < dataPoints[i]["real_pos"] + 50)
                                    dataPoints[i].selected = true;
                            }
                        }
                    }

                    else {
                        for (var i in dataPoints)
                            dataPoints[i].selected = false;
                    }
                    }
                drawPoints(dataPoints, gdl);
                }

        })
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);
}



//-------------------------------------------------------------
//          DATALINE 2
//-------------------------------------------------------------
// partSeq - portion of sequence's abundance to draw
// numNucleotides - number of nucleotides to draw (?) (should be equal to partSeq.length?)
// point - starting point of the drawing
// sizePattern - size of the searched pattern (to highlight over the whole seq)
function dataLine_2(seq, numNucleotides, point, sizePattern, dataName)
{
    var globalDL2=getLevel(dataName, "DL2")
    var globalSeq=getLevel(dataName, "Seq");
    if(globalDL2.cv.nameSVG!="lineSeq2_0")//modify a little the margins if it's not the first one
        globalDL2.cv.margin.bottom=8;

    // DRAWING DATALINE 2
    //-------------------------------------------------
    var startSeq = globalDL2.startSeq;
    var endSeq   = globalDL2.endSeq;

    var gap=(((endSeq-startSeq)/globalSeq.ws-sizePattern)/2)*globalSeq.ws;

    globalDL2.scale=1;
    if(Math.floor((sizePattern*globalSeq.ws)/(endSeq-startSeq)) > 1)
        globalDL2.scale= Math.ceil((sizePattern*globalSeq.ws)/(endSeq-startSeq));


    // We use the core function
    dataLine_core(globalDL2,globalSeq, seq, 1/globalDL2.scale, 0, numNucleotides, startSeq,  point, sizePattern);


    // DRAWING NUCLEOTIDES (DATALINE 3)
    //-------------------------------------------------

    globalDL2.cv.svg.selectAll("."+globalDL2.cv.classSVG+".line").on("mouseover", function (d)
    {
        // DRAWING BRACKETS AND LOCATION
        // Removes previous info about location
        d3.selectAll(".dl2.hover.point").remove();
        d3.selectAll(".dl2.hover.vline").remove();
        d3.selectAll(".dl2.hover.label").remove();

         for(i in dl2) {
            var dl2i=dl2[i];
            var dd=dl2i.cv.svg.selectAll("."+dl2i.cv.classSVG+".line")[0][0]["__data__"];

            var line_x0 = d3.event.pageX - dl2i.cv.margin.left;

            var yScale = d3.scale.linear()
                .domain([globalSeqs[i].min, globalSeqs[i].max])
                .range([dl2i.cv.dim.height, 0]);

            var line_y0 = -1;
            var pointXY = null;

            for (var i in dd)        //This loops over the proper data line
             if (dd[i].pos == line_x0) {
                 line_y0 = yScale(dd[i].value);
                 pointXY = dd[i];
                 break;
             }
            // Determine the width of the brackets (note: approximate)
            var widthText = getTextWidth("A", "12px Courier New");
            var len = (dimDL.width / (widthText * 2))/dl2i.cv.scaleServScreen ;
            drawLocation(dl2i, line_x0, line_y0, pointXY);
            /*if(dl2i.cv.nameSVG==globalDL2.cv.nameSVG)
                drawBrackets(globalDL2, line_x0-len, line_x0+len, line_y0);*/
        }
    /*
        // DRAWING NUCLEOTIDES
        //var start = startSeq + line_x0-len; // taking the left bracket as start
        var start = (startSeq + (line_x0 -len)*dl2i.cv.scaleServScreen ) ; // taking the left bracket as start
        Server.nucleotides(drawNucleotides, globalSeq.track, startSeq, endSeq, start, point, globalSeq.dataName);
   */ });

    if(globalSeqs[0].dataName==dataName)//only on the first line
    {
        //DRAW SHIFT AND ZOOM
        drawShiftZoom(globalDL2);

        var positions={}
        positions[globalSeq.track]=[]
        positions[globalSeq.track].push(point)

        // DRAWING ANNOTATIONS
        Server.annotationsPoint(drawAnnotations, positions, "[\"any\"]", globalDL2.cv.dim.width*2, "center", "False", globalSeq.dataName);
    }
    // We confirm that we have finished
    globalDL2.drawn = true;

}

function getSelectedDataPoint()
    {
    var dataPoint=null;
    var i=0
    //for(var i in dl1)
        {
        if (dl1[i].cv.svg.selectAll("." + dl1[i].cv.classSVG + ".point.pressed").data().length > 0)
            {//Take the first (and only) selected point
            dataPoint = dl1[i].cv.svg.selectAll("." + dl1[i].cv.classSVG + ".point.pressed").data()[0];
            return dataPoint;
            }
        }
    return dataPoint;
    }

function drawAllDataLine_2(seqs, start, end)
    {
    var dataPoint=getSelectedDataPoint();
    for(var i in seqs)
        {
        var globalDL2=getLevel(i,"DL2");
        var globalSeq=getLevel(i,"Seq");
        globalDL2.startSeq = start;
        globalDL2.endSeq = end;
        //dataLine_2(seqs[i],seqs[i].partSeq.length,dataPoint.real_pos,dataPoint.size,globalSeq.dataName)
            dataLine_2(seqs[i],seqs[i].partSeq.length,dataPoint.real_pos,     globalDL2.sizePattern    ,globalSeq.dataName)
        }
    }


/**
 * Draw an annotation gene line below the first DL2 lane.
 * @param annotations
 */
function drawAnnotations(annotations)
{
    var globalDL2=dl2[0];

    var startSeq = globalDL2.startSeq;
    var endSeq = globalDL2.endSeq;
    var factor = globalDL2.cv.scaleSeqScreen;
    var lineHeight = globalDL2.cv.dimAnnotation.lineHeight;

    //gene line
    globalDL2.cv.svg.append("g")
        .selectAll(".dl2.annotation")
        .data(annotations)
        .enter().append("rect")
        .attr('class', 'dl2 annotation')
        .attr("transform", "translate(0," + (globalDL2.cv.dim.height+20) + ")")
        .attr("x", function(d)
        {
        return Math.max(0,(d.s-startSeq)/factor)
        })
        .attr("y", function(d){  var y=0; y=d.ss=="+"?0:15; y+=(d.t=="gene")?5:0; return y;})
        .attr("width", function(d)
            {
            return Math.max(0,Math.min(globalDL2.cv.dim.width-Math.max(0,(d.s-startSeq)/factor)+1, Math.max(0,(d.e-Math.max(startSeq, d.s))/factor)+1))
            })
        .attr("fill", function(d) {
            if (d.t.indexOf("RNA") >= 0)    return "#ffccff";
            if (d.t.indexOf("UTR") < 0)    return "#bbffbb";
            return "#ddffdd";
            })

        .attr("height", function(d){
            var h=0;
            if (d.t=="gene")// || d.t.indexOf("UTR")>=0)
                h=lineHeight*.25;
            else if(d.t=="ORF" || d.t=="CDS" || d.t.indexOf("_gene")>=0)
                h = lineHeight;
            return h});

    //gene arrow
    globalDL2.cv.svg.append("g")
        .selectAll(".dl2.annotation.arrow")
        .data(annotations)
        .enter().append("polygon")
        .attr('class', 'dl2 annotation arrow')
        .attr("transform", "translate(0," + (globalDL2.cv.dim.height+20) + ")")
        .attr("fill", function(d) {
            if (d.t.indexOf("RNA") >= 0)    return "#ffccff";
            return "#ddffdd";
        })

        .attr("points", function(d){
            if(d.t.indexOf("gene")<0 || d.t=="ORF")
                return "";

            var y0=d.ss=="+"?0:15;
            var x0=Math.max(0,(d.s-startSeq)/factor);
            var path="";
            if(d.ss=="+")
                x0+=Math.min(globalDL2.cv.dim.width-Math.max(0,(d.s-startSeq)/factor), Math.max(0,(d.e-Math.max(startSeq, d.s))/factor));
            if(d.ss=="+")
                path=x0+","+y0+ " " +x0+", "+(y0+lineHeight)+" "+(x0+lineHeight *.5)+","+(y0+lineHeight *.5);
            else
                path=x0+","+y0+ " " +x0+", "+(y0+lineHeight)+" "+(x0-lineHeight *.5)+","+(y0+lineHeight *.5);
            if(d.e<startSeq || d.s>startSeq+globalDL2.cv.dim.width)
                return "";
            return path;
        });

    //draw crackers if the annotation goes beyond screen limits
    globalDL2.cv.svg.append("g")
        .selectAll(".dl2.annotation.cracker")
        .data(annotations)
        .enter().append("path")
        .attr('class', 'dl2 annotation cracker')
        .attr("transform", "translate(0," + (globalDL2.cv.dim.height+20) + ")")
        .attr("d", function(d){
            var cw=3;//cracker width and height
            var ch=10;
            var x0=0; var y0=1;
            if (d.s < startSeq)
            {
                x0+=5;
                if(d.ss=="-")
                    y0+=15;
                return "M"+(x0+cw *.5)+","+y0+"L"+x0+","+(y0+ch*.5)+"L"+(x0+cw)+","+(y0+ch*.5)+"L"+x0+","+(y0+ch);
            }
            if (d.e > endSeq)
            {
                if (d.ss == "-")
                    y0 += 15;
                x0+=globalDL2.cv.dim.width-10;
                return "M"+(x0+cw *.5)+","+y0+"L"+x0+","+(y0+ch*.5)+"L"+(x0+cw)+","+(y0+ch*.5)+"L"+x0+","+(y0+ch);
            }
            return "";
        });

    var genes=[];
    for (var key in annotations)
        {
        if (annotations[key]["t"] == "gene" || annotations[key]["t"] == "ORF")
            {
            var g=annotations[key]["id"];
            genes.push(g);
            }
        }


    Server.Gene2GO(setGOterms,genes);


    // Mouseover tip and show the information of genes
    var globalDL1=dl1[0];
    var tip = d3.tip()
        .attr('class', globalDL2.cv.classSVG+' annotation tip')
        .html(function(d,i)
        {
            if(d["text"]==undefined) {
                var text = "<b>id: </b>" + d["id"]+"<br><b>position: </b>" + d["s"] + "-" + d["e"] + " (" + d["ss"] + ")"+
                    "<br><b>type: </b>"+d["t"];

                if (globalDL1.annotations[d["id"]] != undefined) {
                    if (globalDL1.annotations[d["id"]]["goterms"] != undefined) {
                        text += "<br><b>GO terms:</b>";
                        for (var term in globalDL1.annotations[d["id"]]["goterms"]) {
                            var goterm = globalDL1.annotations[d["id"]]["goterms"][term];
                            if (goterm != undefined & goterm.length > 40)    goterm = goterm.substring(0, 40) + "...";
                            text += "<br>  " + goterm;
                            }
                    }
                }
                d["text"]=text;
            }

            return d["text"];
        })
        .offset(function(d){
            return [0,25]}) //above
        ;

    // Calls tip
    globalDL2.cv.svg.call(tip);

    //gene labels
    globalDL2.cv.svg.append("g")
        .selectAll(".dl2.annotation.label")
        .data(annotations)
        .enter().append("text")
        .attr('class', 'dl2 annotation label')
        .attr("transform", "translate(0," + (globalDL2.cv.dim.height+20) + ")")
        .attr("x", function(d) {
            return Math.max(10,(d.s-startSeq+6)/factor);
        })
        .attr("y", function(d){   return d.ss=="+"?10:25})
        .text(function(d){
            var res="";
            if(d.e>startSeq)
                if(d.t.indexOf("gene")>-1 || d.t=="ORF")
                    res=d.n;
            return res;
        })
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide)
        .on('click', function(d){
            tip.hide();

            var startSeq = d.s-50;
            var endSeq   = d.e+50;
            var numNucleotidesDraw=endSeq-startSeq;

            var dataPoint= getSelectedDataPoint();

            Server.getAllPartSeq(drawAllDataLine_2, startSeq, endSeq, dataPoint, dl1,dl2,globalSeqs, 0, {});

        });

}



//-------------------------------------------------------------
//          NUCLEOTIDES (DATALINE 3)
//-------------------------------------------------------------

function removeNucleotides()
    {
        var nameSVG     = globalDL3.cv.nameSVG;
        var svg=d3.select("#"+nameSVG)
        if(svg!=undefined)
            svg.remove();
    }
/**
 *
 * @param start starting position from the mouse hovering
 * @param point point in the sequences retrieved corresponding to our visualization
 * @param nuc all the nucleotides visualized (not only the ones matching the pattern but all the ones in lane2)
 */
function drawNucleotides(start, point, nuc)
{
    // Get information of dataLine
    var nameSVG     = globalDL3.cv.nameSVG;
    var classSVG    = globalDL3.cv.classSVG;
    var width       = globalDL3.cv.dim.width;
    var height      = globalDL3.cv.dim.height;
    var margin      = globalDL3.cv.margin;
    var GDL3Seqs    = globalDL3.sequences;
    var letterWidth = getTextWidth("A", "12px Courier New");

    if(GDL3Seqs.seqs==undefined)    //Not drawn in case that there's only 1 seq or sequences with different size in pattern (by now)
        return;

    // First, we delete the image, if this exist
    var image = $("#"+nameSVG);
    if (image.length) { image.empty(); }
    else
    {
        var core=document.getElementById("core");
        image=document.createElement("div");
        image.setAttribute("id", nameSVG);
        core.appendChild(image);
    }


    // Image SVG: image
    var svg = d3.select("#"+nameSVG)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


    // The image SVG: scale text
    svg.append("g")
        .attr("class", classSVG+" scale")
        .append("text")
        .text(Math.round(letterWidth) + " : 1")
        .attr("x", width-margin.left)
        .attr("y", margin.top-5);

    var seq=undefined
    var motloc=undefined
    if(GDL3Seqs.seqs!=undefined)
        {
        // Draw context brackets
        var x0 = margin.left - 50;
        drawBrackets(globalDL3, x0 - 5, width + 5, globalDL3.cv.bracketHeight * .5);

        seq=GDL3Seqs.seqs[point];      //searched seq (do not confuse with wholeSeq)
        motloc=GDL3Seqs.locations[point]; //TODO: point has not the real position!! --> point.pos+conversion
        }



    var globalDL2=dl2[0];

    var startWholeSeq = globalDL2.startSeq;
    var wholeSeq = nuc.seq;
    start=Math.round(start);
    startWholeSeq=Math.round(startWholeSeq);
    wholeSeq=wholeSeq.toUpperCase();
    var displace=Math.round(Math.max(0,start-point)); //if left bracket is below the first nucleotide of the seq we set to 0
    var numLetters=Math.floor(dimDL.width/letterWidth);

    //Draw the sequence
    var letters = []
    for (var i = start-startWholeSeq; i < start-startWholeSeq+numLetters; i++)
        letters.push(wholeSeq[i]);


    svg.selectAll("mainNucleotides")
        .data(letters)
        .enter()
        .append("text")
        .text(function (d) {
            return d
        })
        .attr("class", classSVG + " text")
        .attr("fill", function (d, i) {
            if (i+start>=point && i+start<point+seq.length)
                return "darkred";
            else
                return "steelblue";
        })
        .style("font-weight", function (d, i) {
            if (start+i>=point+motloc && start+i<point+motloc+GDL3Seqs.motifs[point].length)
                return "bold";
            else
                return "";
        })

        .attr("x", function (d, i) {
            return (x0 + letterWidth * i)
        })
        .attr("y", marginDL.top * 1.7);


    //Draw the remaining sequences (if available):
    if(GDL3Seqs!=undefined) {
        var cont = 1;
        var letterHeight = 12;
        var separator = 5;
        var shownKeys = [];
        for (var key in GDL3Seqs.seqs) {
            if (key != point) {
                shownKeys.push(key);
                var letters0 = []
                var j = 0
                for (var i = 0; i < wholeSeq.length; i++)
                    if (i < point - startWholeSeq || i >= point - startWholeSeq + seq.length)
                        letters0.push("-");
                    else letters0.push(GDL3Seqs.seqs[key][j++]);

                var letters = []
                for (var i = start - startWholeSeq; i < start - startWholeSeq + numLetters; i++)
                    letters.push(letters0[i]);

                var motloc = GDL3Seqs.locations[key];

                svg.selectAll("Nucleotides_" + key)
                    .data(letters)
                    .enter()
                    .append("text")
                    .text(function (d) {
                        return d
                    })
                    .attr("class", classSVG + " text")
                    .attr("fill", function (d, i) {
                        if (i + start >= point && i + start < point + seq.length)
                            return "darkred";
                        else
                            return "steelblue";
                    })
                    .style("font-weight", function (d, i) {
                        if (start + i >= point + motloc && start + i < point + motloc + GDL3Seqs.motifs[point].length)
                            return "bold";
                        else
                            return "";
                    })
                    .attr("x", function (d, i) {
                        return (x0 + letterWidth * i)
                    })
                    .attr("y", marginDL.top * 1.7 + separator + cont * letterHeight);

                cont += 1;
            }
            if (cont > 8)
                break;
        }

        shownKeys.push("consensus");
        //Draw locations
        svg.selectAll("Nucleotides_Text" + key)
            .data(shownKeys)
            .enter()
            .append("text")
            .text(function (d) {
                return d;
            })
            .attr("class", classSVG + " axis y text")
            .attr("x", x0 - 45)
            .attr("y", function (d, i) {
                return marginDL.top * 1.7 + separator * 3.5 + i * letterHeight
            });

        //draw consensus motif
        var letters = []
        for (var i in GDL3Seqs.motifConsensus)
            letters.push(GDL3Seqs.motifConsensus[i])

        var colscale = d3.scale.linear().domain([0, 1]).range(["white", "black"]);

        svg.selectAll("consensus")
            .data(letters)
            .enter()
            .append("text")
            .text(function (d) {
                return d
            })
            .attr("class", classSVG + " consensus")
            .attr("fill", function (d, i) {
                return colscale(GDL3Seqs.motifProfile[d][i])
            })
            .style("font-weight", function (d, i) {
                //if (globalDL3.motifProfile[d][i] > 0.9)
                //    return "bold";
                //else
                return "bold";
            })
            .style("text-decoration", function (d, i) {
                if (GDL3Seqs.motifProfile[d][i] > 0.9)
                    return "underline";
                else
                    return "";
            })
            .attr("x", function (d, i) {
                return (x0 + dimDL.width * .5 - letterWidth * .5 + letterWidth * i)
            })
            .attr("y", marginDL.top * 1.7 + separator * 2 + cont * letterHeight);
        }
    // Save information of dataLine3
    GDL3Seqs.svg = svg;
}


function dataLine_lineOnly(globalDL, globalSeq, seqServ, scaleSeqServ, startSeq, endSeq, initialPoint, point, sizePattern)
    {
    var classSVG    = globalDL.cv.classSVG;
    var nameSVG     = globalDL.cv.nameSVG;
    var height      = globalDL.cv.dim.height;
    var width       = globalDL.cv.dim.width;

    var svg=dl1[0].cv.svg;//Draw all in the same one
    //var svg=globalDL1.cv.svg;//Draw all in the same one



    if(seqServ.minimum!=undefined) //In order to draw variation shade
        {
        var dataMin=scalePoints(seqServ.minimum, startSeq,endSeq, scaleSeqServ, globalDL, globalSeq);
        var dataMax=scalePoints(seqServ.maximum, startSeq,endSeq, scaleSeqServ, globalDL, globalSeq);
        var data=scalePoints(seqServ.partSeq, startSeq,endSeq, scaleSeqServ, globalDL, globalSeq);
        }
    else
        var data=scalePoints(seqServ, startSeq,endSeq, scaleSeqServ, globalDL, globalSeq);

    var ymin = globalSeq.min;
    var ymax=globalSeq.max;

        // Scaling of the axes
    var xScale = d3.scale.linear()
        .domain(d3.extent(data, function (d) { return d.pos; }))  // xmin, xmax
        .range([0, width]);
    var yScale = d3.scale.linear()
        .domain([ymin, ymax])
        .range([height, 0]);

    // Function to draw the line
    var line = d3.svg.line()
        .x(function(d) { return xScale(d.pos); })
        .y(function(d) { return yScale(d.value); });


    // The image SVG: line
    svg.append("g")
        .attr("class", classSVG+" line")
        .datum(data)
        .append("path")
        .attr("stroke",globalDL.cv.color)
        .attr("d", line);

    /*if(globalSeq.correlations==null)
        Server.stats(setStats,globalSeq.dataName);
    else*/
      drawColor(globalDL, globalSeq, svg);

    }

//-------------------------------------------------------------
//          DATALINE CORE AND OTHER FUNCTIONS
//-------------------------------------------------------------
/*
startSeq: first position we want to draw
endSeq: last position we want to draw
seqServ: sequence of points in the above range that we dispose of to draw.
*/
function dataLine_core(globalDL, globalSeq, seqServ, scaleSeqServ, startSeq, endSeq, initialPoint, point, sizePattern)
{
    // Default values
    point       || ( point = 0 );
    sizePattern || ( sizePattern = 0 );

    // Get info about sequence
    var mean = globalSeq.mean;
    var stdev = globalSeq.stdev;



    // Get information of dataLine
    var nameSVG     = globalDL.cv.nameSVG;
    var classSVG    = globalDL.cv.classSVG;
    var width       = globalDL.cv.dim.width;
    var height      = globalDL.cv.dim.height;
    var margin      = globalDL.cv.margin;


    // First, we empy the div, if this exist, or create it
    var image = $("#"+nameSVG);
    if ( image.length) { image.empty(); }
    else
        {
        var core=document.getElementById("core");
        image=document.createElement("div");
        image.setAttribute("id", nameSVG);
        core.appendChild(image);
        }


    /*SCALE: There are three ranges:
        sizeData - Range of nucleotides i want to show
        sizeSeq - Number of data i have in the range of nucleotides to show
        width - Size of the screen

        Based on these ranges, we can define 2/3 scales:

        nucleotides per data -->  seqServ: how many data the server returns me from the whole number of nucleotides
        data per pixel --> servScreen: how many data can i put on each pixel
        nucleotides per pixel --> seqServ: how many actual pixels would i have per pixel

        E.g. the server returns 400K values for a chromosome of 5.5M nucleotides. I have a screen of 1200 pixels.
        Therefore I have: ~13 nucleotides per data (the server returns the average abundance for each 13 nucs)
                          ~355 data would correspond to a single pixel
                          ~355*13=4611 nucleotides would correspond to a single pixel

       Another example: for a small range related to a gene i have 4378 data corresponding to the same number of genes,
       they must be shown again in 1200 pixesl:
       1 nucleotide per data
       3.6 data and nucleotides per pixel

     */
    // Define the scales with that we will work
    var sizeSeq = endSeq-startSeq;
    var startData = Math.floor(startSeq/scaleSeqServ);
    var endData = Math.floor(endSeq/scaleSeqServ);
    var sizeData = endData-startData;

    var scaleServScreen = sizeData/width;
    var scaleSeqScreen=1;
    scaleSeqScreen=scaleSeqServ*scaleServScreen;

    var ymin = globalSeq.min;
    var ymax=globalSeq.max;

    if(seqServ.minimum!=undefined) //In order to draw variation shade
        {
        var dataMin=scalePoints(seqServ.minimum, startSeq,endSeq, scaleSeqServ, globalDL, globalSeq);
        var dataMax=scalePoints(seqServ.maximum, startSeq,endSeq, scaleSeqServ, globalDL, globalSeq);
        var data=scalePoints(seqServ.partSeq, startSeq,endSeq, scaleSeqServ, globalDL, globalSeq);
        }
    else
        var data=scalePoints(seqServ, startSeq,endSeq, scaleSeqServ, globalDL, globalSeq);

    //Scale points for representation

    // Scaling of the axes
    var xScale = d3.scale.linear()
        .domain(d3.extent(data, function (d) { return d.pos; }))  // xmin, xmax
        .range([0, width]);
    var yScale = d3.scale.linear()
        .domain([ymin, ymax])
        .range([height, 0]);

    // Axis labels
    var xAxis = d3.svg.axis().scale(xScale).orient("bottom")
     .tickFormat(function(d) {return roundTickFormat(d,Math.max(scaleSeqScreen, scaleServScreen),initialPoint, startSeq,endSeq)})
     .tickValues(getTicks(sizeSeq, Math.max(scaleSeqScreen, scaleServScreen)));
    var yAxis = d3.svg.axis().scale(yScale).ticks(8).orient("left");

    // Function to draw the line
    var line = d3.svg.line()
        .x(function(d) { return xScale(d.pos); })
        .y(function(d) { return yScale(d.value); });


    // Image SVG: image
    var svg = d3.select("#"+nameSVG)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    globalDL.cv.svg=svg;
    drawGrid();

    // Image SVG: axis x
    svg.append("g")
        .attr("class", classSVG+" x axis")
        .call(xAxis)
        .attr("transform", "translate(0," + height + ")");

    // Image SVG: axis y
    svg.append("g")
        .attr("class", classSVG+" y axis")
        .call(yAxis)
        .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end");

    // Image SVG: data name
    if(globalDL.cv.nameSVG.indexOf("2_")==-1)
        {//DL1 will show just a summary
        /*if(globalSeq.correlations==null)
            Server.stats(setStats,globalSeq.dataName);
        else*/
            drawColor(globalDL, globalSeq, globalDL.cv.svg);
        }
    else
        {
        svg.append("g")
            .attr("class", classSVG + " title")
            .append("text")
            .attr("fill", globalDL.cv.color)//x and y should be in a translate clause if rotation is to be applied
            .attr("transform", "translate(-30," + (globalDL.cv.dim.height * .5 + getTextWidth(globalSeq.dataName, "13px sans-serif") * .5) + ") rotate(-90)")
            .text(globalSeq.dataName);
        }

    //Shadow line in case seqServ provides that info
    if(seqServ.minimum!=undefined) //In order to draw variation shade
        {
            var dataShadow=dataMax;

            for(var i=dataMin.length-1;i>=0;i--)
                {
                var p={};
                p.value = dataMin[i].value;
                p.pos=i;
                dataShadow.push(p);
                }


            svg.append("g")
                .attr("class", classSVG+" lineshadow")
                .datum(dataShadow)
                .append("path")
                .attr("d", line)
                .attr("fill", globalDL.cv.color) //NUEVO
                .attr("stroke", "none");
        }

    // The image SVG: line
    svg.append("g")
        .attr("class", classSVG+" line")
        .datum(data)
        .append("path")
        .attr("stroke",globalDL.cv.color)
        .attr("d", line);





    // The image SVG: highlighted line
    if(point != 0) {
        var startHighlighted = Math.floor(((point - initialPoint) / scaleSeqServ) / scaleServScreen);
        var endHighlighted = Math.floor(((sizePattern * globalSeq.ws) / scaleSeqServ) / scaleServScreen);

        if (startHighlighted + endHighlighted > 0)
            {
            var data2 = data.slice(Math.max(0,startHighlighted), startHighlighted + endHighlighted);

            // The image SVG: line
            svg.append("g")
                .attr("class", classSVG + " line hl")
                .datum(data2)
                .append("path")
                .attr("d", line);
            }
    }


    // The image SVG: scale text
    svg.append("g")
        .attr("class", classSVG+" scale")
        .append("text")
        .text("1 : "+Math.max(Math.round(scaleSeqScreen), Math.round(scaleServScreen)))
        .attr("fill","grey")
        .attr("x", width-margin.left)
        .attr("y", margin.top*0.75);


    // Save information of dataLine
    globalDL.cv.scaleSeqScreen = scaleSeqScreen;
    globalDL.cv.scaleServScreen = scaleServScreen;
    globalDL.cv.data = data;
    globalDL.cv.xScale = xScale;
    globalDL.cv.yScale = yScale;
}


function drawColor(globalDL, globalSeq, svg)
    {
    //var number=parseInt(globalDL.cv.nameSVG.replace(/^.*_/,""))
    var number=globalSeq.id;
    /*var cor=0;
    var cont=0;
    for(var i in globalSeq.correlations)
        {
        cor+=parseFloat(globalSeq.correlations[i]);
        cont+=1
        }
    cor/=cont;
    cor=cor.toFixed(2);*/

    svg.append("g")
        .append("rect")
        .attr("class", globalDL.cv.classSVG +" rect")
        .attr("x",-45)
        .attr("y",10+number*12)
        .attr("fill", globalDL.cv.color)
        .append("svg:title")
        //.text(globalSeq.dataName+"\n  "+globalSeq.numReplicates+" replicates\n  "+cor+" mean correlation");
        .text(globalSeq.dataName);
    }
/**
 * Draws shift (< >) and zoom (+-) elements on a dataline2, plus their logic
 * @param globalDL2
 */
function drawShiftZoom(globalDL2)
    {
    //shift R
    globalDL2.cv.svg.append("g")
        .attr('class', 'dl2 shift')
        .append("text")
        .attr("x", globalDL2.cv.dim.width-5)
        .attr("y", globalDL2.cv.dim.height*.5)
        .text(">")
        .on('mouseover', function(){d3.select(this).attr('fill', '#999')} )
        .on('mouseout', function(){d3.select(this).attr('fill', '#eee')} )
        .on("click",function(){shift(200)});

    //shift L
    globalDL2.cv.svg.append("g")
            .attr('class', 'dl2 shift')
            .append("text")
            .attr("x", 5)
            .attr("y", globalDL2.cv.dim.height*.5)
            .text("<")
            .on('mouseover', function(){d3.select(this).attr('fill', '#999')} )
            .on('mouseout', function(){d3.select(this).attr('fill', '#eee')} )
            .on("click",function(){shift(-200)});


    //zoom + (2x)
    globalDL2.cv.svg.append("g")
        .attr('class', 'dl2 zoom')
        .append("text")
        .attr("x", globalDL2.cv.dim.width-65)
        .attr("y", 8)
        .text("-")
        .on("click",function(){currentZoom=currentZoom*2; zoomLocal(2)});


        //zoom - (0.5x)
    globalDL2.cv.svg.append("g")
        .attr('class', 'dl2 zoom')
        .append("text")
        .attr("x", globalDL2.cv.dim.width-20)
        .attr("y", 8)
        .text("+")
        .on("click",function(){currentZoom=currentZoom*.5; zoomLocal(.5)});

    }

function zoom(amount, start0, end0)
    {
    var start;
    var end;
    var interval=(end0-start0);
    var newInterval=interval*amount;
    var globalDL2=getLevel(globalSeqs[0].dataName, "DL2")
    var scale=Math.max(Math.round(globalDL2.cv.scaleSeqScreen), Math.round(globalDL2.cv.scaleServScreen));


    if(amount>=1)//zoom out
    {
        start = start0 - (newInterval - interval) * .5;
        end = end0 + (newInterval - interval) * .5;
    }
    else{
        if(scale>1) //zoom in up to scale 1:1
        {
            start = start0 + (interval - newInterval) * .5;
            end = end0 - (interval - newInterval) * .5;
        }
        else//zooming in after max zoom in performed, just return
            return;
    }
    var dataPoint= getSelectedDataPoint();
        Server.getAllPartSeq(drawAllDataLine_2, Math.round(start), Math.round(end), dataPoint, dl1,dl2,globalSeqs, 0, {});
    }
//Zooms current view in DL2 by the factor in amount.
function zoomLocal(amount)
    {
    var globalDL2=getLevel(globalSeqs[0].dataName, "DL2")
    zoom(amount, globalDL2.startSeq, globalDL2.endSeq);
    }

//Shifts current view in DL2 by the number of nucleotides in amount.
//Calls an iterative sever requests which finally calls the drawing data line methods
function shift(amount)
    {
    var globalDL2=getLevel(globalSeqs[0].dataName, "DL2")
    var start=Math.round(globalDL2.startSeq+amount);
    var end=Math.round(globalDL2.endSeq+amount);
    var dataPoint= getSelectedDataPoint();
    Server.getAllPartSeq(drawAllDataLine_2, start, end, dataPoint, dl1,dl2,globalSeqs, 0, {});
    }

function getTicks(sizeSeq, scaleSeqScreen)
{
    var ticks = [];

    var numZeroes=0;

    while((sizeSeq/Math.pow(10,numZeroes))>10)
    {
        numZeroes += 1
    }

    if(sizeSeq/Math.pow(10,numZeroes)<4) numZeroes-=1;

    var factorLabel=Math.pow(10, numZeroes);

    //To avoid tick cluttering when multiples of 10 are too many (>20)
    while(sizeSeq/factorLabel>20)
        factorLabel+=factorLabel;

    for(var i=0; factorLabel <= (sizeSeq-i*factorLabel) ;i++)
    {
        ticks.push((factorLabel+i*factorLabel)/scaleSeqScreen);
    }

    return ticks;
}

function scalePoints(seqServ, startSeq, endSeq, scaleSeqServ, globalDL,globalSeq)    // Create the array
    {
    // Define the scales with that we will work
    var sizeSeq = endSeq-startSeq;
    var startData = Math.floor(startSeq/scaleSeqServ);
    var endData = Math.floor(endSeq/scaleSeqServ);
    var sizeData = endData-startData;

    var scaleServScreen = sizeData/globalDL.cv.dim.width;
    var scaleSeqScreen=1;
    scaleSeqScreen=scaleSeqServ*scaleServScreen;

    var data = [];

    for (var i = startData, k = 0; i < endData-Math.ceil(scaleServScreen); i = i + scaleServScreen, k++) {
        var average = 0;
        var numValues = 0;
        for (var j = i; j < i + scaleServScreen; j++)
            {
            if(seqServ[Math.round(j)]!=NaN) {
                average += seqServ[Math.round(j)];
                numValues++;
                }
            else
                console.log("NaN Here!!!")
            }
        if (numValues != 0) average = average / numValues;

        if (average < globalSeq.min) average = globalSeq.min;
        if (average > globalSeq.max) average = globalSeq.max;

        data.push({pos: k, value: average})
        }
    return data;
    }

function roundTickFormat(d, scaleSeqScreen, initialPoint, start, end)
{
    var tickLabel = initialPoint+d*scaleSeqScreen;

    if (tickLabel != "0")
    {
        if(tickLabel/1000000 >= 1 && (end-start)>10000000)//higher than 10M ranges with 1M labels
            tickLabel=Math.round(tickLabel/1000000)+"M";
        else if(tickLabel/1000>=1 && (end-start)>10000) //higher than 10K ranges with 1K labels
            tickLabel=Math.round(tickLabel/1000)+"K";
    }

    return tickLabel;
}


function getTextWidth(text, font)
{
    // Re-use canvas object for better performance
    var canvas = getTextWidth.canvas || (getTextWidth.canvas = document.createElement("canvas"));
    var context = canvas.getContext("2d");
    context.font = font;
    var metrics = context.measureText(text);
    return metrics.width;
}

/**
 * For DL2 hovering
 * @param globalDL
 * @param x0
 * @param y0
 */
function drawLocation(globalDL, x0, y0, point)
{

    // Get information of dataLine and coordinates
    var nameSVG = globalDL.cv.nameSVG;
    var classSVG = globalDL.cv.classSVG;
    var height = globalDL.cv.bracketHeight;
    var desp;

    // POINT
    d3.select("#"+nameSVG)
        .select("svg")
        .select("g")
        .append("circle")
        .attr("cx",x0)
        .attr("cy",y0)
        .attr("class", classSVG+" hover point");

    // LABEL (ABUNDANCE) (SHADOW)
    d3.select("#"+nameSVG)
        .select("svg")
        .select("g")
        .append("text")
        .attr("x",x0-getTextWidth(""+point.value.toFixed(2), "10px sans-serif")*.5)
        .attr("y",y0-5)
        .text(point.value.toFixed(2))
        .attr("class", classSVG+" hover label shadow");

    // LABEL (ABUNDANCE)
    d3.select("#"+nameSVG)
        .select("svg")
        .select("g")
        .append("text")
        .attr("x",x0-getTextWidth(""+point.value.toFixed(2), "10px sans-serif")*.5)
        .attr("y",y0-5)
        .text(point.value.toFixed(2))
        .attr("class", classSVG+" hover label");

    //LABEL (POSITION) (SHADOW)
    var label=d3.select("#"+nameSVG)
        .select("svg")
        .select("g")
        .append("text")
        .attr("x",x0-getTextWidth(""+(globalDL.startSeq+point.pos), "10px sans-serif")*.5)
        .attr("y",globalDL.cv.dim.height+15)
        .text(globalDL.startSeq+point.pos)
        .attr("class", classSVG+" hover label shadow")

    // LABEL (POSITION)
    var label=d3.select("#"+nameSVG)
        .select("svg")
        .select("g")
        .append("text")
        .attr("x",x0-getTextWidth(""+(globalDL.startSeq+point.pos), "10px sans-serif")*.5)
        .attr("y",globalDL.cv.dim.height+15)
        .text(globalDL.startSeq+point.pos)
        .attr("class", classSVG+" hover label")

    return;
}

function drawBrackets(globalDL, left_x0, right_x0, y0)
    {
    // Get information of dataLine and coordinates
    var nameSVG = globalDL.cv.nameSVG;
    var classSVG = globalDL.cv.classSVG;
    var height = globalDL.cv.bracketHeight;
    var x0;
    var desp;
    var bracket;

    // Removes brackets
    d3.selectAll("."+classSVG+".bracket")
        .remove();


    // LEFT BRACKET
    x0   = left_x0;
    desp = 5;
    bracket="M "+(x0+desp)+", "+(y0-height/2)+
        " L "+(x0)+", "+(y0-height/2)+
        " L "+(x0)+", "+(y0+height/2)+
        " L "+(x0+desp)+", "+(y0+height/2);

    d3.select("#"+nameSVG)
        .select("svg")
        .select("g")
        .append("path")
        .attr("d", bracket)
        .attr("class", classSVG+" bracket");


    // RIGHT BRACKET
    x0   = right_x0;
    desp = -5;
    bracket="M "+(x0+desp)+", "+(y0-height/2)+
        " L "+(x0)+", "+(y0-height/2)+
        " L "+(x0)+", "+(y0+height/2)+
        " L "+(x0+desp)+", "+(y0+height/2);

    d3.select("#"+nameSVG)
        .select("svg")
        .select("g")
        .append("path")
        .attr("d", bracket)
        .attr("class", classSVG+" bracket");
    }

function saveSequences(response)
{
    var sequences = [];

    if(response!=null) {

        sequences.seqs = response.seqs;
        sequences.alignment = response.alignment;
        sequences.aconsensus = response.aconsensus;
        sequences.aprofile = response.aprofile;
        sequences.consensus = response.consensus;
        sequences.profile = response.profile;

        sequences.motifs = response.motifs;
        sequences.locations = response.locations;
        sequences.motifConsensus = response.motifConsensus;
        sequences.motifProfile = response.motifProfile;
    }

    globalDL3.sequences = sequences;
}

function setStats(stats, dataName)
    {
    var dl=getLevel(dataName,"Seq")
    dl.numReplicates=stats["numReplicates"]
    dl.correlations=stats["correlation"]

    var dl1=getLevel(dataName,"DL1", true);
    var svg=getLevel(dataName,"DL1").cv.svg;
    drawColor(dl1,dl,svg);
    }

function setAnnotations(gis, annotations)
{
    var globalDL1=dl1[0];
    globalDL1.gis = gis;  // list of gene ids (NOTE: just a string, might contain duplicates, by now unused)
    if (!$.isEmptyObject(annotations))
    {
        globalDL1.annotations = annotations;
    }
}


function setGOterms(genes, goterms)
    {
    var globalDL1=dl1[0];
    for(var i in genes)
        {
        var gene=genes[i];
        if (globalDL1.annotations[gene] == undefined)
            globalDL1.annotations[gene] = {};
        globalDL1.annotations[gene]["goterms"] = goterms[gene];
        }
    }

/** DRAW GRID
 * For datalines 1 and 2 in case the option to see how bands are determined is selected
 */
function drawGrid() {

    if ($("#paramGrid").is(':checked') == false)    //No order to draw
    {
        for(var i in dl1)
            {
            var globalDL1 = dl1[i];
            if (globalDL1.drawn & globalDL1.cv.svg!=undefined)
                {
                globalDL1.cv.svg.selectAll(".dl1.band").remove();
                globalDL1.cv.svg.selectAll(".dl1.band.letter").remove();
                }
            }
        for(var i in dl2)
            {
            var globalDL2=dl2[i]
            if (globalDL2.drawn & globalDL2.cv.svg!=undefined)
                {
                globalDL2.cv.svg.selectAll(".dl1.band").remove();
                globalDL2.cv.svg.selectAll(".dl1.band.letter").remove();
                }
            }
        return;
    }

    for(var i in dl1)
    {
        var globalDL1 = dl1[i];
        var globalDL2 = dl2[i];
        var globalSeq = globalSeqs[i];

        var bands = globalSeq.bins;
        var y0 = 0;
        var margin = globalDL1.cv.margin;
        //var scale = (globalDL1.cv.dim.graphHeight - globalDL1.cv.margin.top - globalDL1.cv.margin.top) / (Math.max(globalSeq.max, bands[bands.length - 1]) - bands[0])
        var scale = (dimDL.height) / (Math.max(globalSeq.max, bands[bands.length - 1]) - bands[0])

        var bw = []; //band widths
        for (var i = bands.length - 1; i >= 1; i--) {
            var top;
            if (i == bands.length - 1)
                top = globalSeq.max;
            else
                top = bands[i];
            top = bands[i];
            console.log("Adding band width for " + bands[i] + "-" + bands[i - 1] + " = " + scale * (top - bands[i - 1]));
            bw.push(Math.max(0, scale * (top - bands[i - 1])));
        }

        //1) BANDS IN LINE 1
        if (globalDL1.cv.svg!= undefined && globalDL1.drawn && globalDL1.cv.svg.selectAll(".dl1.band")[0].length <= 1)    //If not already drawn
        {
            globalDL1.cv.svg.selectAll("bands")
                .data(bw)
                .enter()
                .append("rect")
                .attr('class', 'dl1 band')
                .attr("x", 0)
                .attr("y", function (d, i) {
                    if (i < bands.length - 1) {
                        y0 += d;
                        return y0 - d;
                    }
                    else
                        return 0;
                })
                .attr("fill", function (d, i) {
                    if (i % 2 == 0)
                        return "#ffc477";
                    else
                        return "white";
                })
                .attr("width", dimDL.graphWidth - margin.right - margin.left)
                .attr("height", function (d) {
                    return d;
                });

            y0 = 0;
            var code = 97 + bw.length;

            // Draw band letters
            globalDL1.cv.svg.selectAll("band_letters")
                .data(bw)
                .enter()
                .append("text")
                .attr('class', 'dl1 band letter')
                .attr('x', 10)
                .attr("y", function (d, i) {
                    y0 += d;
                    return y0 - d + 10;
                })
                .text(function (d, i) {
                    code -= 1;
                    return String.fromCharCode(code);
                });
        }

        //2) BANDS IN LINE 2
        if (globalDL2.cv.svg!= undefined && globalDL2.drawn && globalDL2.cv.svg.selectAll(".dl1.band")[0].length <= 1)    //If not already drawn
        //if (globalDL2.drawn)
        {
            y0 = 0;
            globalDL2.cv.svg.selectAll("bands")
                .data(bw)
                .enter()
                .append("rect")
                .attr('class', 'dl1 band')
                .attr("x", 0)
                .attr("y", function (d, i) {
                    if (i < bands.length - 1) {
                        y0 += d;
                        return y0 - d;
                    }
                    else
                        return 0;
                })
                .attr("fill", function (d, i) {
                    if (i % 2 == 0)
                        //return "#ffc477";
                        return "grey";
                    else
                        return "white";
                })
                .attr("width", dimDL.graphWidth - margin.right - margin.left)
                .attr("height", function (d) {
                    return d;
                });

            y0 = 0;
            var code = 97 + bw.length;

            // Draw occurrences label
            globalDL2.cv.svg.selectAll("band_letters")
                .data(bw)
                .enter()
                .append("text")
                .attr('class', 'dl1 band letter')
                .attr('x', 10)
                .attr("y", function (d, i) {
                    y0 += d;
                    return y0 - d + 10;
                })
                .text(function (d, i) {
                    code -= 1;
                    return String.fromCharCode(code);
                });
        }
    }//for each data line name
}