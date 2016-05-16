/*
 ┌────────────────────────────────────────────────────────────┐
 │ draw.js                                                    │
 ├────────────────────────────────────────────────────────────┤
 │ Description:                                               │
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
    graphHeight : 200,

    width : screen.width-marginDL.right-marginDL.left,
    height : 200-marginDL.top-marginDL.bottom
};

var dimDLAnnotation =
{
    x0:5,
    y0:7,
    lineHeight : 12
};


var globalSeq =
{
    seqServ : null,
    fullLength : null,
    mean : null,
    stdev : null,
    ws : null,
    scaleSeqServ : null,
    tracks : null
};


// Chromosome level
var globalDL1 =
{
    // Core variables (used in the core)
    cv :
    {
        nameSVG : "lineSeq",
        classSVG : "dl1",

        svg : null,
        scaleSeqScreen : null,
        scaleServScreen : null,
        data : null,
        xScale : null,
        yScale : null,

        margin :
        {
            top: marginDL.top,
            bottom:marginDL.bottom,
            left:marginDL.left,
            right:marginDL.right
        },
        dim :
        {
            graphWidth : dimDL.graphWidth,
            graphHeight : dimDL.graphHeight,
            width : dimDL.width,
            height : dimDL.height
        },
        dimAnnotation :
        {
            x0 : dimDLAnnotation.x0,
            y0 : dimDLAnnotation.y0,
            lineHeight : dimDLAnnotation.lineHeight
        }
    },

    drawn : false,

    seqPoints : null
};


// Nucleosome level
var globalDL2 =
{
    // Core variables (used in the core)
    cv :
    {
        nameSVG : "lineSeq2",
        classSVG : "dl2",

        svg : null,
        scaleSeqScreen : null,
        scaleServScreen : null,
        data : null,
        xScale : null,
        yScale : null,

        margin :
        {
            top: marginDL.top,
            bottom:marginDL.bottom,
            left:marginDL.left,
            right:marginDL.right
        },
        dim :
        {
            graphWidth : dimDL.graphWidth,
            graphHeight : dimDL.graphHeight,
            width : dimDL.width,
            height : dimDL.height
        },
        dimAnnotation :
        {
            lineHeight : dimDLAnnotation.lineHeight
        },
        bracketHeight: 40
    },

    startSeq : null,
    endSeq : null
};


// Nucleotide level
var globalDL3 =
{
    // Core variables (used to draw)
    cv :
    {
        nameSVG : "lineSeq3",
        classSVG : "dl3",

        svg : null,
        data : null,


        margin :
        {
            top: marginDL.top,
            bottom:marginDL.bottom,
            left:marginDL.left,
            right:marginDL.right
        },
        dim :
        {
            graphWidth : dimDL.graphWidth,
            graphHeight : dimDL.graphHeight,
            width : dimDL.width,
            height : dimDL.height
        },
        bracketHeight: dimDL.height
    },

    sequences :
    {
        seqs: null,
        alignment: null,
        aconsensus: null,
        aprofile: null,
        consensus: null,
        profile: null,

        motifs: null,
        locations: null,
        motifConsensus: null,
        motifProfile: null
    }
};



//-------------------------------------------------------------
//          DATALINE 1
//-------------------------------------------------------------


function dataLine_1(track, fullLength, seqServ, startSeq, endSeq, maxSize, mean, stdev, ws)
{
    var startTime = new Date();

    // Save information of the genomic sequence
    globalSeq.ws            = ws;
    globalSeq.track         = track; // chromosome or wig track it refers to
    globalSeq.seqServ       = seqServ;
    globalSeq.mean          = mean;
    globalSeq.stdev         = stdev;
    globalSeq.scaleSeqServ  = 1;
    if(Math.floor(fullLength/maxSize) >= 1)
        globalSeq.scaleSeqServ = Math.floor(fullLength/maxSize);


    // We use the core function
    dataLine_core(false, globalDL1,
        globalSeq.seqServ, globalSeq.scaleSeqServ, startSeq, endSeq, startSeq);


    // We confirm that we have finished
    globalDL1.drawn = true;

    if(DEBUG_GBV) console.log("Time spent dataLine1: "+ (new Date()-startTime)+"ms");
}

function dataLine_1_drawPoints(allPoints, tracks, track, sizePattern, numNucleotidesDraw)
{
    var points=JSON.parse(allPoints[track]);
    if(DEBUG_GBV) console.log("\ndataLine_1_drawPoints(): "+points.length+" points");

    // We calculate the points found in the sequence
    var seqPoints=[];
    for(var i=0; i<points.length;i++)
    {
        seqPoints.push(points[i] * globalSeq.ws);
    }

    // Now, with seqPoints, we calculate the points found in the dataLine
    var dataPoints=[];
    for(i=0; i<seqPoints.length;i++)
    {
        var dataPoint = Math.round(seqPoints[i]/globalDL1.cv.scaleSeqScreen);
        dataPoints.push({real_pos: seqPoints[i], pos: (globalDL1.cv.data)[dataPoint].pos, value: (globalDL1.cv.data)[dataPoint].value});
    }


    // Mouseover tip and drawing the corresponding line
    var tip = d3.tip()
        .attr('class', globalDL1.cv.classSVG+' point-tip')
        .offset([30, 0]) // [top, left] to center the tip
        .html(function(d,i)
        {
            var point = d.real_pos;

            // We round to the nearest hundred from sizePattern. E.g. 270 --> 300,  540 --> 500
            var focusLine = Math.round(sizePattern*globalSeq.ws/100)*100;

            var startSeq = point-(numNucleotidesDraw/2)+(focusLine/2);
            var endSeq   = point+(numNucleotidesDraw/2)+(focusLine/2);

            // Save information of dataLine2
            globalDL2.startSeq = startSeq;
            globalDL2.endSeq = endSeq;

            // Draw dataLine2
            Server.getPartSeq(dataLine_2, globalSeq.track, startSeq, endSeq, numNucleotidesDraw, point, sizePattern);

            return "<strong>"+d3.format(",")(point) + ":</strong> " + d3.format(".2f")(d.value); // e.g. "307,770 : 0.79"
        });




        // Calls tip
    globalDL1.cv.svg.call(tip);


    // Remove all points (previous) and occurrences label TODO: cambiar a globalDL1.cv sin más?
    globalDL1.cv.svg.selectAll(".point")
        .remove();
    globalDL1.cv.svg.selectAll(".search-label")
        .remove();


    // Create all new points
    globalDL1.cv.svg.selectAll(".data")
        .data(dataPoints)
        .enter()
        .append("circle")
        .attr('class', globalDL1.cv.classSVG+' point')
        .attr('cx', function(d) { return globalDL1.cv.xScale(d.pos); })
        .attr('cy', function(d) { return globalDL1.cv.yScale(d.value); })
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);


    // Keep colored the selected point
    $("."+globalDL1.cv.classSVG+".point").bind( "mouseover", function()
    {
        $("."+globalDL1.cv.classSVG+".point").attr('class', globalDL1.cv.classSVG+' point');
        $(this).attr('class', globalDL1.cv.classSVG+' point pressed');
    });


    // Draw occurrences label
    globalDL1.cv.svg.append("g")
        .attr("class", globalDL1.cv.classSVG+" search-label")
        .append("text")
        .attr('x', 385)
        .attr('y', -2)
        .text("matches:");

     var matches=[]
     for(var i in tracks) {
     var pp = JSON.parse(allPoints[tracks[i]]);
     matches.push(pp.length);
         }
     globalDL1.cv.svg.selectAll(".searchLabels")
     .data(matches)
     .enter()
     .append("text")
     .attr("class", globalDL1.cv.classSVG+" search-label")
     .attr('x', function(d,i){return 440+i*30;})
     .attr('y', -2)
     .style("font-weight", function (d, i) {
         if(tracks[i]==track) return "bold";
         return "";
     })
     .text(function(d){
         return d+""
     });

    // Save information of dataLine1
    globalDL1.seqPoints = seqPoints;


    //GET SEQUENCES, MOTIFS, (ALIGNMENT), CONSENSUS
    //NOTE: alignment takes more than 1s if there's >50 sequences! (using the fastest method: kalign)
    for( var j=0;j<points.length;j++)
        points[j]*=globalSeq.ws;
    Server.nucProfile(saveSequences, globalSeq.track, "["+points+"]", sizePattern*globalSeq.ws);

}

function dataLine_1_drawEnrichment(enrichment)
{
    // Remove all goterm (previous)
    globalDL1.cv.svg.selectAll("."+globalDL1.cv.classSVG+".goterm")
        .remove();


    // Adjust the scale
    var goSize = d3.scale.log().base(10)
        .domain([1,10e-6])  // p-values
        .clamp(true)
        .rangeRound([6, 14]); // min and max letter size

    // We create a array with goterms
    var goterms=[];
    for(var k in enrichment)
    {
        if(enrichment[k].go_name != undefined)
            goterms.push(enrichment[k]);
        else
            console.log("dataLine_1_drawEnrichment(): Error! => GO term "+k+" not found (possibly outdated OBO file?");

    }
    if(DEBUG_GBV) console.log("dataLine_1_drawEnrichment(): "+goterms.length+" enriched terms");


    // Mouseover tip and show the information of goterm
    var tip = d3.tip()
        .attr('class', globalDL1.cv.classSVG+' goterm-tip')
        .offset([45, 0])   // [top, left] to center the tip
        .html(function(d,i)
        {
            var text = "p-value: "+ d3.format(".2e")(d.pval)+"<br>"+
                d.ngis+"/"+ d.ngo +" genes"+"<br>";
            return text;

            /*if(fixed)
             {
             console.log("Fixed!");
             for (var g in d.gis)
             text += d.gis[g] + "<br>";
             }*/
            //TODO: highlight the gis in the track visible now (infrastructure ready on globalDL1.iannotations and enrichment.gis)
            //Such infrastructure requires pre-load of every position annotations: quite expensive in time at the search
            /*console.log("Related positions:")
             for(var g in d.gis)ls

             console.log(d.gis[g]+"\t"+globalDL1.iannotations[globalSeq.track][d.gis[g]]);*/
            //This is asking now to the go terms annotations
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
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);
    //.on('click', function(){tip.fixed?fix.show(true):fix.show(false); tip.fixed=!tip.fixed;});
}


//-------------------------------------------------------------
//          DATALINE 2
//-------------------------------------------------------------

function dataLine_2(partSeq, numNucleotides, point, sizePattern)
{
    // DRAWING DATALINE 2
    //-------------------------------------------------
    var startSeq = globalDL2.startSeq;
    var endSeq   = globalDL2.endSeq;


    // We use the core function
    dataLine_core(false, globalDL2,
        partSeq, 1, 0, numNucleotides, startSeq,
        point, sizePattern);



    // DRAWING ANNOTATIONS
    //-------------------------------------------------
    Server.annotationsGenes(dataLine_2_drawAnnotationLine, point,"[\"any\"]",globalDL2.cv.dim.width, "center", globalSeq.track, "False");


    // DRAWING NUCLEOTIDES (DATALINE 3)
    //-------------------------------------------------
    globalDL2.cv.svg.selectAll("."+globalDL2.cv.classSVG+".line").on("mouseover", function (d)
    {

        // DRAWING BRACKETS
        var line_x0 = d3.event.layerX-globalDL2.cv.margin.left;
        var line_y0 = d3.event.layerY-dimDL.height*2;

        // Determine the width of the brackets (note: approximate)
        var widthText = getTextWidth("A", "12px Courier New");
        var len = dimDL.width/(widthText*2);

        drawBrackets(globalDL2, line_x0-len, line_x0+len, line_y0);

        // DRAWING NUCLEOTIDES
        var start = startSeq+line_x0-len; //taking the left bracket as start
        //var start=startSeq+d3.event.layerX-marginDL.left; //by now just taking the matched sequence.
        Server.nucleotides(drawNucleotides, globalSeq.track, startSeq, endSeq, start, point);//and the nucleotides for the third lane
    });
}

function dataLine_2_drawAnnotationLine(annotations)
{
    if(false) console.log("\ndataLine_2_drawAnnotationLine(): there are annotations");


    var startSeq = globalDL2.startSeq;
    var endSeq = globalDL2.endSeq;
    var factor=globalDL2.cv.scaleSeqScreen*globalDL2.cv.scaleServScreen; // De momento queremos que el factor sea 1:1 siempre
    var lineHeight = globalDL2.cv.dimAnnotation.lineHeight;

    //console.log(annotations);
    //console.log(factor);


    //gene line
    globalDL2.cv.svg.append("g")
        .selectAll(".dl2.annotation")
        .data(annotations)
        .enter().append("rect")
        .attr('class', 'dl2 annotation')
        .attr("transform", "translate(0," + (globalDL2.cv.dim.height+20) + ")")
        .attr("x", function(d)
        {  return Math.max(0,(d.start-startSeq)*factor) })
        .attr("y", function(d){  var y=0; y=d.sense=="+"?0:15; y+=(d.type=="gene" || d.type=="transcript")?5:0; return y;})
        .attr("width", function(d){return Math.max(0,Math.min(globalDL2.cv.dim.width-Math.max(0,(d.start-startSeq)*factor), Math.max(0,d.end-Math.max(startSeq, d.start)*factor)))})
        .attr("height", function(d){
            var h=lineHeight;
            if (d.type=="gene" || d.type=="transcript")
                h=lineHeight*.25;
            return h});

    //gene arrow
    globalDL2.cv.svg.append("g")
        .selectAll(".dl2.annotation.arrow")
        .data(annotations)
        .enter().append("polygon")
        .attr('class', 'dl2 annotation arrow')
        .attr("transform", "translate(0," + (globalDL2.cv.dim.height+20) + ")")
        .attr("points", function(d){
            var y0=d.sense=="+"?0:15;
            var x0=Math.max(0,(d.start-startSeq)*factor);
            if(d.sense=="+")
                x0+=Math.min(globalDL2.cv.dim.width-Math.max(0,(d.start-startSeq)*factor), Math.max(0,d.end-Math.max(startSeq, d.start)*factor));
            if(d.sense=="+")
                path=x0+","+y0+ " " +x0+", "+(y0+lineHeight)+" "+(x0+lineHeight *.5)+","+(y0+lineHeight *.5);
            else
                path=x0+","+y0+ " " +x0+", "+(y0+lineHeight)+" "+(x0-lineHeight *.5)+","+(y0+lineHeight *.5);
            if(d.end<startSeq || d.start>startSeq+globalDL2.cv.dim.width)
                return "";
            if(d.type.indexOf("gene")>-1)
                return path;
            else
                return "";
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
            if (d.start < startSeq)
            {
                x0+=5;
                if(d.sense=="-")
                    y0+=15;
                return "M"+(x0+cw *.5)+","+y0+"L"+x0+","+(y0+ch*.5)+"L"+(x0+cw)+","+(y0+ch*.5)+"L"+x0+","+(y0+ch);
            }
            if (d.end > endSeq)
            {
                if (d.sense == "-")
                    y0 += 15;
                x0+=globalDL2.cv.dim.width-10;
                return "M"+(x0+cw *.5)+","+y0+"L"+x0+","+(y0+ch*.5)+"L"+(x0+cw)+","+(y0+ch*.5)+"L"+x0+","+(y0+ch);
            }
            return "";
        });



    //gene labels
    globalDL2.cv.svg.append("g")
        .selectAll(".dl2.annotation.label")
        .data(annotations)
        .enter().append("text")
        .attr('class', 'dl2 annotation label')
        .attr("transform", "translate(0," + (globalDL2.cv.dim.height+20) + ")")
        .attr("x", function(d) {
            return Math.max(10,(d.start-startSeq+6)*factor);
        })
        .attr("y", function(d){   return d.sense=="+"?10:25})
        .text(function(d){
            res="";
            if(d.end>startSeq)
                if(d.type.indexOf("gene")>-1)
                    res=d.id;
            return res;
        });
}



//-------------------------------------------------------------
//          NUCLEOTIDES (DATALINE 3)
//-------------------------------------------------------------

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


    // First, we delete the image, if this exist
    var image = $("#"+nameSVG);
    if ( image.length) { image.empty(); }


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


    // Draw context brackets
    var x0 = margin.left - 50;
    drawBrackets(globalDL3, x0 - 5, width + 5, height/2);



    var startWholeSeq = globalDL2.startSeq;
    var wholeSeq = nuc.seq;
    start=Math.round(start);
    startWholeSeq=Math.round(startWholeSeq);
    wholeSeq=wholeSeq.toUpperCase();
    var displace=Math.round(Math.max(0,start-point)); //if left bracket is below the first nucleotide of the seq we set to 0
    var numLetters=Math.floor(dimDL.width/letterWidth);
    var seq=GDL3Seqs.seqs[point];      //searched seq (do not confuse with wholeSeq)

    //Draw the sequence
    var letters = []
    for (var i = start-startWholeSeq; i < start-startWholeSeq+numLetters; i++)
        letters.push(wholeSeq[i]);
    var motloc=GDL3Seqs.locations[point]; //TODO: check -1 values in locations

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


    //Draw the remaining sequences:
    var cont = 1;
    var letterHeight = 12;
    var separator = 5;
    var shownKeys=[];
    for (var key in GDL3Seqs.seqs) {
        if (key != point) {
            shownKeys.push(key);
            var letters0 = []
            var j=0
            for (var i = 0; i<wholeSeq.length;i++)
                if(i<point-startWholeSeq || i>=point-startWholeSeq+seq.length)
                    letters0.push("-");
                else letters0.push(GDL3Seqs.seqs[key][j++]);

            var letters=[]
            for (var i = start-startWholeSeq; i < start-startWholeSeq+numLetters; i++)
                letters.push(letters0[i]);

            var motloc=GDL3Seqs.locations[key];

            svg.selectAll("Nucleotides_"+key)
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
                .attr("y", marginDL.top * 1.7 + separator + cont * letterHeight);

            cont += 1;
        }
        if (cont > 8)
            break;
    }

    shownKeys.push("consensus");
    //Draw locations
    svg.selectAll("Nucleotides_Text"+key)
        .data(shownKeys)
        .enter()
        .append("text")
        .text(function (d) { return d; })
        .attr("class", classSVG+" axis y text")
        .attr("x", x0-45)
        .attr("y", function(d,i){ return marginDL.top * 1.7 + separator*3.5 + i * letterHeight});

    //draw consensus motif
    var letters = []
    for (var i in GDL3Seqs.motifConsensus)
        letters.push(GDL3Seqs.motifConsensus[i])

    var colscale=d3.scale.linear().domain([0,1]).range(["white", "black"]);

    console.log("scale done");
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
        .style("text-decoration", function(d,i){
            if (GDL3Seqs.motifProfile[d][i] > 0.9)
                return "underline";
            else
                return "";})
        .attr("x", function (d, i) {
            return (x0 +dimDL.width *.5-letterWidth *.5 + letterWidth * i)
        })
        .attr("y", marginDL.top * 1.7 + separator * 2 + cont * letterHeight);


    // Save information of dataLine3
    GDL3Seqs.svg = svg;
}





//-------------------------------------------------------------
//          DATALINE CORE AND OTHER FUNCTIONS
//-------------------------------------------------------------

function dataLine_core(DEBUG, globalDL,
                       seqServ, scaleSeqServ, startSeq, endSeq, initialPoint,
                       point, sizePattern)
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


    // First, we delete the image, if this exist
    var image = $("#"+nameSVG);
    if ( image.length) { image.empty(); }


    // Define the scales with that we will work
    var sizeSeq = endSeq-startSeq;
    var startData = Math.floor(startSeq/scaleSeqServ);
    var endData = Math.floor(endSeq/scaleSeqServ);
    var sizeData = endData-startData;

    var scaleServScreen = sizeData/width;
    var scaleSeqScreen=1;
    if(sizeSeq>width)   // width/bps compression
    {
        scaleSeqScreen = sizeSeq/width; // i.e. nucleotides per pixel
    }
    if(DEBUG) console.log("dataLine("+classSVG+"): startSeq: "+startSeq+" - endSeq: "+endSeq+" - sizeSeq: "+sizeSeq+" - width(pixels): "+width);
    if(DEBUG) console.log("            startData: "+startData+" - endData: "+endData+" - sizeData: "+sizeData);
    if(DEBUG) console.log("            scaleSeqScreen (nucleotides/pixel): "+scaleSeqScreen+" - scaleServScreen: "+scaleServScreen+" - scaleSeqServ: "+scaleSeqServ);


    // Calculate 'y minimum' and 'y maximum'
    var ymin = 0;
    var ymax = mean+(3*stdev);
    if(mean-(3*stdev) > 0)
        ymin = mean-3*stdev;


    // Create the array
    var data = [];
    for(var i=startData, k=0; i < endData; i=i+scaleServScreen, k++)
    {
        var average = 0;
        var numValues = 0;
        for (var j=i; j < i+scaleServScreen; j++)
            if (seqServ[Math.round(j)] >= 0)
            {
                average += seqServ[Math.round(j)];
                numValues++;
            }
        if(numValues != 0) average = average/numValues;

        if (average < ymin) average = ymin;
        if (average > ymax) average = ymax;

        data.push({pos: k, value: average})
    }


    // Scaling of the axes
    var xScale = d3.scale.linear()
        .domain(d3.extent(data, function (d) { return d.pos; }))  // xmin, xmax
        .range([0, width]);
    var yScale = d3.scale.linear()
        .domain([ymin, ymax])
        .range([height, 0]);

    // Axis labels
    var xAxis = d3.svg.axis().scale(xScale).orient("bottom")
        .tickFormat(function(d) {return roundTickFormat(d,scaleSeqScreen,initialPoint, startSeq,endSeq)})
        .tickValues(getTicks(sizeSeq, scaleSeqScreen));
    var yAxis = d3.svg.axis().scale(yScale).orient("left");

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

    // The image SVG: line
    svg.append("g")
        .attr("class", classSVG+" line")
        .datum(data)
        .append("path")
        .attr("d", line);


    // The image SVG: highlighted line
    if(point!=0)
    {
        var startHighlighted = Math.floor(((point-initialPoint)/scaleSeqServ)/scaleServScreen);
        var endHighlighted = Math.floor(((sizePattern*globalSeq.ws)/scaleSeqServ)/scaleServScreen);

        var data2 = data.slice(startHighlighted, startHighlighted+endHighlighted);

        // The image SVG: line
        svg.append("g")
            .attr("class", classSVG + " line hl")
            .datum(data2)
            .append("path")
            .attr("d", line);
    }

    // The image SVG: scale text
    svg.append("g")
        .attr("class", classSVG+" scale")
        .append("text")
        .text("1 : "+Math.round(scaleSeqScreen))
        .attr("x", width-margin.left)
        .attr("y", margin.top*0.75);


    // Save information of dataLine
    globalDL.cv.svg = svg;
    globalDL.cv.scaleSeqScreen = scaleSeqScreen;
    globalDL.cv.scaleServScreen = scaleServScreen;
    globalDL.cv.data = data;
    globalDL.cv.xScale = xScale;
    globalDL.cv.yScale = yScale;
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

    for(var i=0; factorLabel <= (sizeSeq-i*factorLabel) ;i++)
    {
        ticks.push((factorLabel+i*factorLabel)/scaleSeqScreen);
    }

    return ticks;
}

function roundTickFormat(d, scaleSeqScreen, initialPoint, start, end)
{
    var tickLabel = initialPoint+d*scaleSeqScreen;

    if (tickLabel != "0")
    {
        if(tickLabel/1000000 >= 1 && (end-start)>10000000)
            tickLabel=Math.round(tickLabel/1000000)+"M";
        else if(tickLabel/1000>=1 && (end-start)>10000)
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
        .attr("class", classSVG+" bracket")


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
        .attr("class", classSVG+" bracket")
}

function saveSequences(response)
{
    var sequences = [];

    sequences.seqs=response.seqs;
    sequences.alignment=response.alignment;
    sequences.aconsensus=response.aconsensus;
    sequences.aprofile=response.aprofile;
    sequences.consensus=response.consensus;
    sequences.profile=response.profile;

    sequences.motifs=response.motifs;
    sequences.locations=response.locations;
    sequences.motifConsensus=response.motifConsensus;
    sequences.motifProfile=response.motifProfile;

    globalDL3.sequences = sequences;
}

/*
function setAnnotations(gis, annotations)
{
    globalDL1.gis=gis;                 //list of gene ids
    if(annotations!=null)
    {
        globalDL1.annotations = annotations;//position -> element id (element can be a gene or other entity)
        globalDL1.iannotations = {}       //gene id -> position
        for (var p in annotations) {
            globalDL1.iannotations[p] = {}
            for (var e in annotations[p]) {
                for (var i = 0; i < annotations[p][e].length; i++)
                    if (annotations[p][e][i]["type"] == "gene")
                        globalDL1.iannotations[p][annotations[p][e][i]["id"]] = e
            }
        }
    }
}
*/




