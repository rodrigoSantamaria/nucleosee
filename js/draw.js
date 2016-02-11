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
        y0:5
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





function dataLine1(seqServ, startSeq, endSeq, fullLength, maxSize, mean, stdev, ws)
{
    // Get information of dataLine1
    var nameSVG=globalDL1.nameSVG;
    var classSVG=globalDL1.classSVG;
    var width=globalDL1.dim.graphWidth;
    var height=globalDL1.dim.graphHeight;

    // Save information of sequence
    globalSeq.seqServ=seqServ;
    globalSeq.mean=mean;
    globalSeq.stdev=stdev;
    globalSeq.ws=ws;
    globalSeq.scaleSeqServ = 1;
    if(Math.floor(fullLength/maxSize) >= 1)
        globalSeq.scaleSeqServ = Math.floor(fullLength/maxSize);


//    var result = drawDataLine(nameSVG, classSVG, width, height,marginDL,
    var result = drawDataLine(nameSVG, classSVG,
        seqServ, globalSeq.scaleSeqServ, startSeq, startSeq, endSeq);


    // Save information of dataLine1
    //globalDL1.width = result.width;
    //globalDL1.height = result.height;
    globalDL1.scaleSeqScreen = result.scaleSeqScreen;
    globalDL1.scaleServScreen = result.scaleServScreen;
    globalDL1.data = result.data;
    globalDL1.x = result.x;
    globalDL1.y = result.y;
    globalDL1.svg = result.svg;

    // We confirm that we have finished
    globalDL1.drawn = true;
}

function drawEnrichment(enrichment)
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

        var tip = d3.tip()
            .attr('class', 'dl1 goterm-tip')
            .offset([45, 0])
            //.style("opacity",.8)
            .attr("fill","blue")
            .html(function(d,i)
            {
                return "p-value: "+ d3.format(".2e")(d.pval) +"<br>" + d.ngis+"/"+ d.ngo +" genes<br>";
            });

        globalDL1.svg.call(tip);


        var dx=dimAnnotation.x0;

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
                    dx+=width;//+dimAnnotation.x0;
                    return dx-width; }
                )
            .attr('y', function(d) { return dimAnnotation.y0 })
            .attr('font-size', function(d) { //console.log("Tamaño de "+ d.go_name+": "+d.pval+" "+goSize(d.pval)+"px");
                return goSize(d.pval)+"px" })
            .text(function(d){
                //console.log(d.go_name);
                return " "+d.go_name+" ·"})
            .on('mouseover', tip.show)
            .on('mouseout', tip.hide);

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
        dataPoints.push({pos: (globalDL1.data)[dataPoint].pos, value: (globalDL1.data)[dataPoint].value});
    }


    // Mouseover tip and drawing the corresponding line
    var tip = d3.tip()
        .attr('class', 'dl1 point-tip')
        .offset([30, 0])
        .html(function(d,i)
        {
            dataLine2(seqPoints[i], sizePattern, numNucleotides);
            return "<strong>"+d3.format(",")(seqPoints[i]) + ":</strong> " + d3.format(".2f")(d.value);
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
        .on('mouseover', tip.show)
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
        .attr('y', 20)
        .text(function(d){return d+" occurrences"})

    // Save information of dataLine1
    globalDL1.seqPoints = seqPoints;
}


function dataLine2(point, sizePattern, numNucleotides)
{
    // Get information of dataLine2
    var nameSVG=globalDL2.nameSVG;
    var classSVG=globalDL2.classSVG;
    var width=globalDL2.dim.graphWidth;
    var height=globalDL2.dim.graphHeight;


    // We round to the nearest hundred from sizePattern
    var focusLine = Math.round(sizePattern*globalSeq.ws/100)*100;

    var startSeq = point-(numNucleotides/2)+(focusLine/2);
    var endSeq   = point+(numNucleotides/2)+(focusLine/2);

    // Get the part of sequence that we need to draw
    var result = Server.getPartSeq(startSeq,endSeq);


   // var dataLine = drawDataLine(nameSVG, classSVG, width, height, globalDL2.paddingX, globalDL2.paddingY,
   // var dataLine = drawDataLine(nameSVG, classSVG, width, height, marginDL,
    var dataLine = drawDataLine(nameSVG, classSVG,
        result.partSeq, 1, startSeq, 0, numNucleotides, // scaleSeqServ is 1 because is a part of sequence
                                point, sizePattern);

    //RODRIGO
    var annotations = Server.annotationsGenes("["+point+"]","[\"any\"]",globalDL2.dim.width, "center");
    console.log("annotations result hola "+annotations);
    if(annotations.hasOwnProperty(point))
        var annotLine = drawAnnotationLine(dataLine, annotations[point], startSeq, endSeq);
    //RODRIGO


    // Save information of dataLine2
    //globalDL2.width = dataLine.width;
    //globalDL2.height = dataLine.height;
    globalDL2.scaleSeqScreen = dataLine.scaleSeqScreen;
    globalDL2.scaleServScreen = dataLine.scaleServScreen;
    globalDL2.data = dataLine.data;
    globalDL2.x = dataLine.x;
    globalDL2.y = dataLine.y;
    globalDL2.svg = dataLine.svg;
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

        //TODO: far below/above cracker?

        //gene labels
    dataLine.svg.append("g")
        .selectAll(".dl2.annotation.label")
        .data(annotations)
        .enter().append("text")
        .attr('class', 'dl2 annotation label')
        .attr("transform", "translate(0," + (globalDL2.dim.height+20) + ")")
        .attr("x", function(d) {
            /*if(d.sense=="-")
               return Math.min(globalDL2.dim.width-20,(-startSeq+ d.end+3)*factor);
            else
                return Math.max(0,((d.start-startSeq)-100)*factor);*/
            return Math.max(0,(d.start-startSeq+6)*factor);
            })
        .attr("y", function(d){   return d.sense=="+"?10:25})
        .text(function(d){
            //console.log(d.id+" "+ d.type +" ["+ d.start+", "+ d.end+"]"+" "+ d.sense);
            res=""
            if(d.end>startSeq)
                if(d.type.indexOf("gene")>-1)
                    res=d.id;
            return res
        });
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
    if(DEBUG_GBV) console.log("\ndataLine(): startSeq: "+startSeq+" - endSeq: "+endSeq+" - sizeSeq: "+sizeSeq+" - graphWidth: "+dimDL.width);
    if(DEBUG_GBV) console.log("            startData: "+startData+" - endData: "+endData+" - sizeData: "+sizeData);
    if(DEBUG_GBV) console.log("            scaleSeqScreen (nucleotides/pixel): "+scaleSeqScreen+" - scaleServScreen: "+scaleServScreen+" - scaleSeqServ: "+scaleSeqServ);


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
