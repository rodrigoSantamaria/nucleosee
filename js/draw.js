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
    graphHeight : 200,
    graphWidth : screen.width,

    width : screen.width-marginDL.right-marginDL.left,
    height : 200-marginDL.top-marginDL.bottom
};

var dimAnnotation =
    {
        height : 12,
        x0:5,
        y0:7
    }

var globalSeq =
{
    seqServ : null,
    fullLength : null,
    mean : null,
    stdev : null,
    ws : null,
    scaleSeqServ : null
};


//chromosome level
var globalDL1 =
{
    dim : dimDL,
    margin : marginDL,

    drawn : false,

    nameSVG : "lineSeq",
    classSVG : "dl1",

    scaleSeqScreen : null,
    scaleServScreen : null,
    data : null,
    x : null,
    y : null,
    svg : null,
    seqPoints : null
};

//nucleosome level
var globalDL2 =
{
    dim : dimDL,
    margin : marginDL,

    nameSVG : "lineSeq2",
    classSVG : "dl2",

    scaleSeqScreen : null,
    scaleServScreen : null,
    data : null,
    x : null,
    y : null,
    svg : null
};

//nucleotide level
var globalDL3 =
{
    dim : dimDL,
    margin : marginDL,

    nameSVG : "lineSeq3",
    classSVG : "dl3",

    scaleSeqScreen : null,
    scaleServScreen : null,
    data : null,
    x : null,
    y : null,
    svg : null,

    alignment: null,
    aconsensus: null,
    aprofile: null,
    seqs: null,
    consensus: null,
    profile: null
};





function dataLine1(seqServ, startSeq, endSeq, fullLength, maxSize, mean, stdev, ws, track)
{
    // Get information of dataLine1
    var nameSVG=globalDL1.nameSVG;
    var classSVG=globalDL1.classSVG;
    var width=globalDL1.dim.graphWidth;
    var height=globalDL1.dim.graphHeight;

    // Save information of the genomic sequence
    globalSeq.seqServ=seqServ;
    globalSeq.mean=mean;
    globalSeq.stdev=stdev;
    globalSeq.ws=ws;
    globalSeq.track=track;//chromosome or wig track it refers to
    globalSeq.scaleSeqServ = 1;
    if(Math.floor(fullLength/maxSize) >= 1)
        globalSeq.scaleSeqServ = Math.floor(fullLength/maxSize);

    var result = drawDataLine(nameSVG, classSVG,
        seqServ, globalSeq.scaleSeqServ, startSeq, startSeq, endSeq);


    // Save information of dataLine1
    globalDL1.scaleSeqScreen = result.scaleSeqScreen;
    globalDL1.scaleServScreen = result.scaleServScreen;
    globalDL1.data = result.data;
    globalDL1.x = result.x;
    globalDL1.y = result.y;
    globalDL1.svg = result.svg;

    // We confirm that we have finished
    globalDL1.drawn = true;
}


function drawEnrichment(enrichment, annotations)
    {
        globalDL1.svg.selectAll(".dl1.goterm")
            .remove();

        var goSize = d3.scale.log().base(10)
            .domain([1,10e-6])  // p-values
            .clamp(true)
            .rangeRound([6, 14]); //letter size

        var goterms=[];
        for(var k in enrichment) {
            if(enrichment[k].go_name!=undefined)
                goterms.push(enrichment[k])
            else
                console.log("Error!: GO term "+k+" not found (possibly outdated OBO file?");

        }

        console.log(goterms.length+" enriched terms");

        var tip = d3.tip()
            .attr('class', 'dl1 goterm-tip')
            .offset([45, 0])
            .html(function(d,i)
            {
                var text="p-value: "+ d3.format(".2e")(d.pval) +"<br>" + d.ngis+"/"+ d.ngo +" genes<br>";
                /*if(fixed)
                    {
                    console.log("Fixed!");
                    for (var g in d.gis)
                        text += d.gis[g] + "<br>";
                    }*/
                //TODO: highlight the gis in the track visible now (infrasructure ready on globalDL1.iannotations and enrichment.gis)
                //Such infrastructure requires pre-load of every position annotations: quite expensive in time at the search
                /*console.log("Related positions:")
                for(var g in d.gis)ls

                    console.log(d.gis[g]+"\t"+globalDL1.iannotations[globalSeq.track][d.gis[g]]);*/
                //This is asking now to the go terms annotations
                return text;
            });

        globalDL1.svg.call(tip);

        var dx=dimAnnotation.x0;
        var dy=[];
        var yline=0;

        console.log("Numero de terminos GO: "+goterms.length);
        globalDL1.svg.selectAll(".dl1.goterm")
            .data(goterms)
            .enter()
            .append("text")
            .attr('class', 'dl1 goterm')
            .attr('x', function(d)
                {
                    var canvas = document.createElement('canvas');
                    var ctx = canvas.getContext("2d");
                    ctx.font = goSize(d.pval)+"px sans-serif";
                    var width = ctx.measureText(" "+d.go_name+" ·").width;

                    dx+=width;
                    if(dx>=dimDL.width) {
                        yline++;
                        dx = dimAnnotation.x0 + width;
                        }
                    dy[d.go_name]=dimAnnotation.y0+globalDL1.dim.height+25+10*yline;
                    console.log(d.go_name+" "+ dy[d.go_name]);
                    return dx-width; }
                )
            .attr('y', function(d) {
                return dy[d.go_name];//*Math.max(1,dx/dimDL.width)
            })
            .attr('font-size', function(d) { //console.log("Tamaño de "+ d.go_name+": "+d.pval+" "+goSize(d.pval)+"px");
                return goSize(d.pval)+"px" })
            .text(function(d){
                //console.log(d.go_name);
                return " "+d.go_name+" ·"})
            .on('mouseover', tip.show)
            .on('mouseout', tip.hide)
            //.on('click', function(){tip.fixed?fix.show(true):fix.show(false); tip.fixed=!tip.fixed;});

    }

function drawPoints(points, sizePattern, numNucleotides)
{
    var seqPoints=[];
    for(var i=0; i<points.length;i++)
    {
        seqPoints.push(points[i] * globalSeq.ws);
    }

    if(DEBUG_GBV) console.log("\ndataLine1_drawPoints():");
    if(DEBUG_GBV) console.log(seqPoints);

    var dataPoints=[];
    for(i=0; i<seqPoints.length;i++)
    {
        var dataPoint = Math.floor(seqPoints[i]/globalDL1.scaleSeqScreen);
        dataPoints.push({real_pos: seqPoints[i], pos: (globalDL1.data)[dataPoint].pos, value: (globalDL1.data)[dataPoint].value});
    }


    // Mouseover tip and drawing the corresponding line
    var tip = d3.tip()
       // .attr('class', 'dl1 point-tip')
        .offset([30, 0])
        .html(function(d,i)
        {
            dataLine2(d.real_pos, sizePattern, numNucleotides);
            return "<strong>"+d3.format(",")(d.real_pos) + ":</strong> " + d3.format(".2f")(d.value);
        });



    // Calls tip
    globalDL1.svg.call(tip);

    // Remove all points
    globalDL1.svg.selectAll(".point")
        .remove();

    // Create all new points
    globalDL1.svg.selectAll(".data")
        .data(dataPoints)
        .enter()
        .append("circle")
        .attr('class', 'dl1 point')
        .attr('cx', function(d) { return globalDL1.x(d.pos) })
        .attr('cy', function(d) { return globalDL1.y(d.value) })
    .on('mouseover', function(d){  d.attr('class', 'dl1 point selected'); tip.show} )
    .on('mouseout', tip.hide);

    globalDL1.svg.selectAll(".search_label")
        .remove();

    // Draw occurrences label
    globalDL1.svg.selectAll(".search_label")
        .data([dataPoints.length])
        .enter()
        .append("text")
        .attr('class', 'dl1 search_label')
        .attr('x', 600)
        .attr('y', 0)
        .text(function(d){return d+" matches"})

    // Save information of dataLine1
    globalDL1.seqPoints = seqPoints;
}


function drawAnnotationLine(dataLine, annotations, startSeq, endSeq)
    {
    var factor=dataLine.scaleSeqScreen*dataLine.scaleServScreen;

        //gene line
    dataLine.svg.append("g")
        .selectAll(".dl2.annotation")
        .data(annotations)
        .enter().append("rect")
        .attr('class', 'dl2 annotation')
        .attr("transform", "translate(0," + (globalDL2.dim.height+20) + ")")
        .attr("x", function(d)
            {  return Math.max(0,(d.start-startSeq)*factor) })
        .attr("y", function(d){  var y=0; y=d.sense=="+"?0:15; y+=(d.type=="gene" || d.type=="transcript")?5:0; return y;})
        .attr("width", function(d){return Math.max(0,Math.min(globalDL2.dim.width-Math.max(0,(d.start-startSeq)*factor), Math.max(0,d.end-Math.max(startSeq, d.start)*factor)))})
        .attr("height", function(d){
            var h=dimAnnotation.height;
            if (d.type=="gene" || d.type=="transcript")
                h=dimAnnotation.height*.25;
            return h});

        dataLine.svg.append("g")
            .selectAll(".dl2.annotation.arrow")
            .data(annotations)
            .enter().append("polygon")
            .attr('class', 'dl2 annotation arrow')
            .attr("transform", "translate(0," + (globalDL2.dim.height+20) + ")")
            .attr("points", function(d){
                //TODO
                var y0=d.sense=="+"?0:15;
                var x0=Math.max(0,(d.start-startSeq)*factor);
                if(d.sense=="+")
                    x0+=Math.min(globalDL2.dim.width-Math.max(0,(d.start-startSeq)*factor), Math.max(0,d.end-Math.max(startSeq, d.start)*factor));
                if(d.sense=="+")
                    path=x0+","+y0+ " " +x0+", "+(y0+dimAnnotation.height)+" "+(x0+dimAnnotation.height *.5)+","+(y0+dimAnnotation.height *.5);
                else
                    path=x0+","+y0+ " " +x0+", "+(y0+dimAnnotation.height)+" "+(x0-dimAnnotation.height *.5)+","+(y0+dimAnnotation.height *.5);
                if(d.end<startSeq || d.start>startSeq+globalDL2.dim.width)
                    return "";
                if(d.type.indexOf("gene")>-1)
                    return path;
                else
                    return "";
            });


    //draw crackers if the annotation goes beyond screen limits
    dataLine.svg.append("g")
        .selectAll(".dl2.annotation.cracker")
        .data(annotations)
        .enter().append("path")
        .attr('class', 'dl2 annotation cracker')
        .attr("transform", "translate(0," + (globalDL2.dim.height+20) + ")")
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
                x0+=globalDL2.dim.width-10;
                return "M"+(x0+cw *.5)+","+y0+"L"+x0+","+(y0+ch*.5)+"L"+(x0+cw)+","+(y0+ch*.5)+"L"+x0+","+(y0+ch);
                }
            return "";
        });



        

        //gene labels
    dataLine.svg.append("g")
        .selectAll(".dl2.annotation.label")
        .data(annotations)
        .enter().append("text")
        .attr('class', 'dl2 annotation label')
        .attr("transform", "translate(0," + (globalDL2.dim.height+20) + ")")
        .attr("x", function(d) {
            return Math.max(10,(d.start-startSeq+6)*factor);
            })
        .attr("y", function(d){   return d.sense=="+"?10:25})
        .text(function(d){
            res=""
            if(d.end>startSeq)
                if(d.type.indexOf("gene")>-1)
                    res=d.id;
            return res
        });
    }


function dataLine2(point, sizePattern, numNucleotides) {
    // Get information of dataLine2
    var nameSVG = globalDL2.nameSVG;
    var classSVG = globalDL2.classSVG;
    var width = globalDL2.dim.graphWidth;
    var height = globalDL2.dim.graphHeight;

    //remove the nucleotide line too
    var imageSVG = $("#lineSeq3");
    if ( imageSVG.length) { imageSVG.empty(); }

    // We round to the nearest hundred from sizePattern
    var focusLine = Math.round(sizePattern * globalSeq.ws / 100) * 100;

    var startSeq = point - (numNucleotides / 2) + (focusLine / 2);
    var endSeq = point + (numNucleotides / 2) + (focusLine / 2);

    // Get the part of sequence that we need to draw
    var result = Server.getPartSeq(startSeq, endSeq, globalSeq.track);//the abundance
    var nuc=Server.nucleotides(startSeq, endSeq, globalSeq.track);//and the nucleotides for the third lane

    // var dataLine = drawDataLine(nameSVG, classSVG, width, height, globalDL2.paddingX, globalDL2.paddingY,
    // var dataLine = drawDataLine(nameSVG, classSVG, width, height, marginDL,
    var dataLine = drawDataLine(nameSVG, classSVG,
        result.partSeq, 1, startSeq, 0, numNucleotides, // scaleSeqServ is 1 because is a part of sequence
        point, sizePattern);

    // <<<  RODRIGO
    // SEQUENCE
    dataLine.svg.selectAll(".dl2.line").on("mouseover", function (d) {

        //DRAWING BRACKETS
        var bh=40;//bracket height
        var x0=d3.event.layerX-marginDL.left;
        var y0=d3.event.layerY-dimDL.height*2+30;
        //console.log(x0+","+y0+" start: "+start);

        //determine the width of the brackets (note: approximate)
        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext("2d");
        ctx.font = "12px Courier New";
        var width = ctx.measureText("A").width;
        var len=dimDL.width/(width*2)

        //var start=startSeq+d3.event.layerX-marginDL.left; //by now just taking the matched sequence.
        var start=startSeq+d3.event.layerX-marginDL.left-len; //taking the left bracket as start

        //replace brackets
        d3.selectAll(".dl2.line.bracket").remove();//removes everything, not only brackets
        drawBracket(nameSVG,"dl2 line",x0-len, y0, true, bh);
        drawBracket(nameSVG,"dl2 line",x0+len, y0, false, bh);

        //DRAWING NUCLEOTIDES
        drawNucleotides("lineSeq3", "dl3", start, point, nuc.seq, startSeq);
    });

    //ANNOTATIONS
    var annotations = Server.annotationsGenes("["+point+"]","[\"any\"]",globalDL2.dim.width, "center", globalSeq.track, "False");
    if(annotations.hasOwnProperty(point))
        var annotLine = drawAnnotationLine(dataLine, annotations[point], startSeq, endSeq);
    //NOTE: this if we compute annotations previously, in the moment of the search (might be time consuming)
    //if(globalDL1.annotations[globalSeq.track].hasOwnProperty(point))
    //    var annotLine = drawAnnotationLine(dataLine, globalDL1.annotations[globalSeq.track][point], startSeq, endSeq);
    // >>>>>> RODRIGO


    // Save information of dataLine2
    globalDL2.scaleSeqScreen = dataLine.scaleSeqScreen;
    globalDL2.scaleServScreen = dataLine.scaleServScreen;
    globalDL2.data = dataLine.data;
    globalDL2.x = dataLine.x;
    globalDL2.y = dataLine.y;
    globalDL2.svg = dataLine.svg;
}

function setSequences(response)
    {
    globalDL3.seqs=response.seqs;
        globalDL3.alignment=response.alignment;
        globalDL3.aconsensus=response.aconsensus;
        globalDL3.aprofile=response.aprofile;
        globalDL3.consensus=response.consensus;
        globalDL3.profile=response.profile;

        globalDL3.motifs=response.motifs;
        globalDL3.locations=response.locations;
        globalDL3.motifConsensus=response.motifConsensus;
        globalDL3.motifProfile=response.motifProfile;
    }

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

/**
 *
 * @param nameSVG Name of the svg tag hosting the image
 * @param classSVG Class of the svg tag hosting the image
 * @param start starting position from the mouse hovering
 * @param point point in the sequences retrieved corresponding to our visualization
 * @param wholeSeq all the nucleotides visualized (not only the ones matching the pattern but all the ones in lane2)
 * @returns {{width: *, height: *, data: *, x0: number, y0: number, svg: *}}
 */
function drawNucleotides(nameSVG, classSVG, start, point, wholeSeq, startWholeSeq)
    {
        start=Math.round(start);
        startWholeSeq=Math.round(startWholeSeq);
        wholeSeq=wholeSeq.toUpperCase();

        var canvas = document.createElement('canvas');
        var ctx = canvas.getContext("2d");
        ctx.font = "12px Courier New";
        var letterWidth = ctx.measureText("A").width;

        var displace=Math.round(Math.max(0,start-point)); //if left bracket is below the first nucleotide of the seq we set to 0

        var numLetters=Math.floor(dimDL.width/letterWidth);
        var seq=globalDL3.seqs[point];      //searched seq (do not confuse with wholeSeq)
        var startSearchedSeq=point

            // First, we delete the image, if this exist
            var imageSVG = $("#" + nameSVG);
            if (imageSVG.length) {
                imageSVG.empty();
            }

            // Image SVG: image
            var svg = d3.select("#" + nameSVG)
                .append("svg")
                .attr("width", dimDL.width + marginDL.left + marginDL.right)
                .attr("height", dimDL.height + marginDL.top + marginDL.bottom)
                .append("g")
                .attr("transform", "translate(" + marginDL.left + "," + marginDL.top + ")");
            //.attr("class", "dl3 ")

            svg.selectAll("text").remove();

            // The image SVG: scale text
            svg.append("text")
                .text(Math.round(letterWidth) + " : 1")
                .attr("class", classSVG + " scale")
                .attr("x", dimDL.width - marginDL.left * 2.5)
                .attr("y", marginDL.top - 5);

            var x0 = marginDL.left - 50;

            //draw context brackets
            drawBracket(nameSVG, classSVG, x0 - 5, dimDL.height * 0.5, true, dimDL.height);
            drawBracket(nameSVG, classSVG, dimDL.width + 5, dimDL.height * 0.5, false, dimDL.height);


            //Draw the sequence
            var letters = []
            for (var i = start-startWholeSeq; i < start-startWholeSeq+numLetters; i++)
                letters.push(wholeSeq[i]);
            var motloc=globalDL3.locations[point]; //TODO: check -1 values in locations

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
                    if (start+i>=point+motloc && start+i<point+motloc+globalDL3.motifs[point].length)
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
            for (var key in globalDL3.seqs) {
                if (key != point) {
                    shownKeys.push(key);
                    var letters0 = []
                    var j=0
                    for (var i = 0; i<wholeSeq.length;i++)
                        if(i<point-startWholeSeq || i>=point-startWholeSeq+seq.length)
                            letters0.push("-");
                        else letters0.push(globalDL3.seqs[key][j++]);

                    var letters=[]
                    for (var i = start-startWholeSeq; i < start-startWholeSeq+numLetters; i++)
                        letters.push(letters0[i]);

                    var motloc=globalDL3.locations[key];

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
                            if (start+i>=point+motloc && start+i<point+motloc+globalDL3.motifs[point].length)
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
            for (var i in globalDL3.motifConsensus)
                letters.push(globalDL3.motifConsensus[i])

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
                    return colscale(globalDL3.motifProfile[d][i])
                })
                .style("font-weight", function (d, i) {
                    //if (globalDL3.motifProfile[d][i] > 0.9)
                    //    return "bold";
                    //else
                        return "bold";
                })
                .style("text-decoration", function(d,i){
                    if (globalDL3.motifProfile[d][i] > 0.9)
                        return "underline";
                    else
                        return "";})
                .attr("x", function (d, i) {
                    return (x0 +dimDL.width *.5-letterWidth *.5 + letterWidth * i)
                })
                .attr("y", marginDL.top * 1.7 + separator * 2 + cont * letterHeight);


        /*    //Draw the consensus sequence (by now not considered useful, maybe with alignments...)
            if (globalDL3.aconsensus == undefined) {
                var letters = []
                for (var i = displace; i < Math.min(displace+numLetters, globalDL3.consensus.length); i++)
                    letters.push(globalDL3.consensus[i])
                svg.selectAll("consensus")
                    .data(letters)
                    .enter()
                    .append("text")
                    .text(function (d) {
                        return d
                    })
                    .attr("class", classSVG + " consensus")
                    .attr("fill", function (d, i) {
                        if (globalDL3.profile[d][displace+i] > 0.5)
                        //if(globalDL3.profile[i]>0.5)
                            return "black";
                        else if (globalDL3.profile[d][displace+i] > 0.25)
                        //else if(globalDL3.profile[i]>0.25)
                            return "lightgrey";
                        else
                            return "white";
                    })
                    .style("font-weight", function (d, i) {
                        if (globalDL3.profile[d][displace+i] > 0.75)
                        //if(globalDL3.profile[i]>0.75)
                            return "bold";
                        else
                            return "";
                    })
                    .attr("x", function (d, i) {
                        return (x0 + letterWidth * i)
                    })
                    .attr("y", marginDL.top * 1.7 + separator * 2 + cont * letterHeight);
            }
            else {
                //draw consensus from the alignment: it's ready but we need to think more about it
                var letters = []
                for (var i = 0; i < globalDL3.aconsensus.length; i++)
                    letters.push(globalDL3.aconsensus[i])
                svg.selectAll("aconsensus")
                    .data(letters)
                    .enter()
                    .append("text")
                    .text(function (d) {
                        return d
                    })
                    .attr("class", classSVG + " consensus")
                    .attr("fill", function (d, i) {
                        if (globalDL3.aprofile[d][i] > 0.5)
                        //if(globalDL3.aprofile[i]>0.5)
                            return "black";
                        else if (globalDL3.aprofile[d][i] > 0.25)
                        //else if(globalDL3.aprofile[i]>0.25)
                            return "lightgrey";
                        else
                            return "white";
                    })
                    .style("font-weight", function (d, i) {
                        if (globalDL3.aprofile[d][i] > 0.75)
                        //if(globalDL3.aprofile[i]>0.75)
                            return "bold";
                        else
                            return "";
                    })
                    .attr("x", function (d, i) {
                        return (x0 + letterWidth * i)
                    })
                    .attr("y", marginDL.top * 1.7 + 15);
            }*/
        globalDL3.svg = svg;
        globalDL3.data=seq;

        // Save information of dataLine1
        return {
            width : dimDL.width,
            height : dimDL.height,
            data : seq,
            x0 : marginDL.left,
            y0 : marginDL.top,
            svg : svg
        }
    }

//function drawDataLine(nameSVG, classSVG, width, height, margin, seqServ, scaleSeqServ, initialPoint, startSeq, endSeq, point, sizePattern)
function drawDataLine(nameSVG, classSVG, seqServ, scaleSeqServ, initialPoint, startSeq, endSeq, point, sizePattern)
    {
    // Default values
    point       || ( point = 0 );
    sizePattern || ( sizePattern = 0 );


    // First, we delete the image, if this exist
    var imageSVG = $("#"+nameSVG);
    if ( imageSVG.length) { imageSVG.empty(); }

    // Get info about sequence
    var mean=globalSeq.mean;
    var stdev=globalSeq.stdev;


    // Define the scales with that we will work
    var startData = Math.floor(startSeq/scaleSeqServ);
    var endData = Math.ceil(endSeq/scaleSeqServ);
    var sizeData = endData-startData;
    var sizeSeq = endSeq-startSeq;


    var scaleServScreen = sizeData/dimDL.width;
    var scaleSeqScreen=1;
    if(sizeSeq>dimDL.width)   // width/bps compression
    {
        scaleSeqScreen=sizeSeq/dimDL.width; // i.e. nucleotides per pixel
    }
    /*if(DEBUG_GBV) console.log("\ndataLine(): startSeq: "+startSeq+" - endSeq: "+endSeq+" - sizeSeq: "+sizeSeq+" - graphWidth: "+dimDL.width);
    if(DEBUG_GBV) console.log("            startData: "+startData+" - endData: "+endData+" - sizeData: "+sizeData);
    if(DEBUG_GBV) console.log("            scaleSeqScreen (nucleotides/pixel): "+scaleSeqScreen+" - scaleServScreen: "+scaleServScreen+" - scaleSeqServ: "+scaleSeqServ);
    */

    // Calculate 'y minimum' and 'y maximum'
    var ymin=0;
    var ymax=mean+3*stdev;
    if(mean-3*stdev>0)
        ymin=mean-3*stdev;


    // Create the array
    var data=[];
    for(var i=startData,k=0;i<endData;i=i+scaleServScreen,k++) {
        var average = 0;
        var numValues = 0;
        for (var j = i; j < i + scaleServScreen; j++)
            if (seqServ[Math.round(j)] >= 0)
            {
                average += seqServ[Math.round(j)];
                numValues++;
            }
        if(numValues != 0) average=average/numValues;

        if (average<ymin) average=ymin;
        if (average>ymax) average=ymax;

        data.push({pos: k, value: average})
    }


    // Scaling of the axes
    var x = d3.scale.linear()
        .domain(d3.extent(data, function (d) { return d.pos; }))  // xmin, xmax
        .range([0, dimDL.width]);
    var y = d3.scale.linear()
        .domain([ymin, ymax])
        .range([dimDL.height, 0]);

    // Axis labels
    var xAxis = d3.svg.axis().scale(x).orient("bottom")
        .tickFormat(function(d) {return roundTickFormat(d,scaleSeqScreen,initialPoint, startSeq,endSeq)})
        .tickValues(getTicks(sizeSeq, scaleSeqScreen));
    var yAxis = d3.svg.axis().scale(y).orient("left");

    // Function to draw the line
    var line = d3.svg.line()
        .x(function(d) { return x(d.pos); })
        .y(function(d) { return y(d.value); });


    // Image SVG: image
    var svg = d3.select("#"+nameSVG)
        .append("svg")
        .attr("width", dimDL.width + marginDL.left + marginDL.right)
        .attr("height", dimDL.height + marginDL.top + marginDL.bottom)
        .append("g")
        .attr("transform", "translate(" + marginDL.left + "," + marginDL.top + ")");

        // Image SVG: axis x
    svg.append("g")
        .attr("class", classSVG+" x axis")
        .call(xAxis)
        .attr("transform", "translate(0," + dimDL.height + ")");

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
        svg.append("path")
            .datum(data2)
            .attr("class", classSVG + " line hl")
            .attr("d", line);
    }

    // The image SVG: scale text
    svg.append("text")
        .text("1 : "+Math.round(scaleSeqScreen))
        .attr("class", classSVG+" scale")
        .attr("x", dimDL.width-marginDL.left*2.5)
        .attr("y", marginDL.top *.75);


    // Save information of dataLine1
    return {
        width : dimDL.width,
        height : dimDL.height,
        scaleSeqScreen : scaleSeqScreen,
        scaleServScreen : scaleServScreen,
        data : data,
        x : x,
        y : y,
        x0 : marginDL.left,  //TODO: this is confusing
        y0 : marginDL.top,
        svg : svg
    }
}

/**
 * Draws a bracket figure
 * @param svg to which we append the brackets
 * @param x0 middle point x of the bracket
 * @param y0 middle point y of the braket vertical line
 * @param left true if it's a left bracket, false if it's right
 * @param height height of the bracket
 */
function drawBracket(nameSVG, classSVG, x0, y0, left, height, clean)
    {
    var desp=left?5:-5;
    var bracket="M "+(x0+desp)+", "+(y0-height *.5)+
        " L "+(x0)+", "+(y0-height *.5)+
        " L "+(x0)+", "+(y0+height *.5)+
        " L "+(x0+desp)+", "+(y0+height *.5);

        if(clean)
            {
                d3.select("#"+nameSVG)
                    .select("svg")
                    .select("g")
                    .selectAll("")
            }
    // The image SVG: line
    d3.select("#"+nameSVG)
        .select("svg")
         .select("g")
        .append("path")
        .attr("d", bracket)
        .attr("class", classSVG+" bracket")
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
