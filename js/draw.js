/*
 ┌────────────────────────────────────────────────────────────┐
 │ draw.js                                                    │
 ├────────────────────────────────────────────────────────────┤
 │ Description:                                               │
 └────────────────────────────────────────────────────────────┘
 */


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
    drawn : false,

    nameSVG : "lineSeq",
    classSVG : "dl1",

    graphHeight : 200,
    graphWidth : screen.width,
    paddingX : 50,
    paddingY : 25,

    width : null,
    height : null,
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
    nameSVG : "lineSeq2",
    classSVG : "dl2",

    graphHeight : 200,
    graphWidth : screen.width,
    paddingX : 50,
    paddingY : 25,

    width : null,
    height : null,
    scaleSeqScreen : null,
    scaleServScreen : null,
    data : null,
    x : null,
    y : null,
    svg : null,
    startSeq : null,
    endSeq : null
};


var globalDL3 =
{
    nameSVG : "lineSeq3",
    classSVG : "dl3",

    graphHeight : 110,  // recommendation: number divisible by 3 and 2 times more dd
    graphWidth : screen.width,
    paddingX : 50,
    paddingY : 25,

    annotations : null
};



function dataLine1(seqServ, startSeq, endSeq, fullLength, maxSize, mean, stdev, ws)
{
    // Get information of dataLine1
    var nameSVG=globalDL1.nameSVG;
    var classSVG=globalDL1.classSVG;
    var width=globalDL1.graphWidth;
    var height=globalDL1.graphHeight;
    var x0=globalDL1.paddingX;
    var y0=globalDL1.paddingY;

    // Save information of sequence
    globalSeq.seqServ=seqServ;
    globalSeq.mean=mean;
    globalSeq.stdev=stdev;
    globalSeq.ws=ws;
    globalSeq.scaleSeqServ = 1;
    if(Math.floor(fullLength/maxSize) >= 1)
        globalSeq.scaleSeqServ = Math.floor(fullLength/maxSize);


    var result = drawDataLine(nameSVG, classSVG, width, height, x0, y0,
                              seqServ, globalSeq.scaleSeqServ, startSeq, startSeq, endSeq);


    // Save information of dataLine1
    globalDL1.width = result.width;
    globalDL1.height = result.height;
    globalDL1.scaleSeqScreen = result.scaleSeqScreen;
    globalDL1.scaleServScreen = result.scaleServScreen;
    globalDL1.data = result.data;
    globalDL1.x = result.x;
    globalDL1.y = result.y;
    globalDL1.svg = result.svg;

    // We confirm that we have finished
    globalDL1.drawn = true;
}


function dataLine1_drawPoints(points, sizePattern, numNucleotides)
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
        if(DEBUG_GBV) console.log("        dataPoint: "+dataPoint);
        dataPoints.push({pos: (globalDL1.data)[dataPoint].pos, value: (globalDL1.data)[dataPoint].value});
    }


    // Mouseover tip and drawing the corresponding line
    var tip = d3.tip()
        .attr('class', 'dl1 point-tip')
        .offset([120, 40])
        .html(function(d,i)
        {
            dataLine2(seqPoints[i], sizePattern, numNucleotides);
            dataLine3(seqPoints[i]);
            return "<strong>" + seqPoints[i] + " position</strong><br>" + Math.round(d.value*100)/100 + " value" + "<br>";
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
        .attr('cx', function(d) { return globalDL1.x(d.pos); })
        .attr('cy', function(d) { return globalDL1.y(d.value); })
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);

    // Save information of dataLine1
    globalDL1.seqPoints = seqPoints;
}




function dataLine3(point)
{
    if (globalDL3.annotations)
    {
        var width   = globalDL3.graphWidth;
        var height  = globalDL3.graphHeight;
        var x0      = globalDL3.paddingX;
        var y0      = globalDL3.paddingY;
        var nameSVG = globalDL3.nameSVG;
        var annots  = globalDL3.annotations[point];

        // Sizes...
        var margin = {top: y0, right: x0, bottom: y0, left: x0};
        width = width - margin.left - margin.right;
        height = height - margin.top - margin.bottom;
        var ySenseLine      = Math.round(height/3);
        var yAntisenseLine  = Math.round(height*2/3);
        var xStartLine = 0;
        var xEndLine  = 0;
        var yLine  = 0;


        // Image SVG: image
        var svg = d3.select("#"+nameSVG)
            .append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
            .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

        // Image SVG: rectangle to know the area where I'm drawing
        svg.append("rect")
            .attr("class", "rect dl3")
            .attr("x", 0)
            .attr("y", 0)
            .attr("width", width)
            .attr("height", height);

        console.log("Inicio: "+globalDL2.startSeq+" - final: "+globalDL2.endSeq);



        /*
            Tipos de anotaciones:
            - gene
            - transcript
            - exon
            - CDS
            - five_prime_UTR
            - three_prime_UTR
            - biological_region
         */


        for(var i=0, j=annots.length; i<j;i++)
        {
            if(annots[i].end>=globalDL2.startSeq && globalDL2.endSeq>=annots[i].start)
            {
                if(i<=9)console.log("[ "+i+"] Start: "+annots[i].start+" - end: "+annots[i].end+" - y es un: "+annots[i].type+" ("+annots[i].id+", "+annots[i].sense+")");
                if(i>9)  console.log("["+i+"] Start: "+annots[i].start+" - end: "+annots[i].end+" - y es un: "+annots[i].type+" ("+annots[i].id+", "+annots[i].sense+")");



                if(annots[i].start-globalDL2.startSeq > 0)
                {
                    xStartLine = Math.round((annots[i].start-globalDL2.startSeq)/globalDL2.scaleSeqScreen);
                }
                else
                {
                    xStartLine = 0;
                }
                if(annots[i].end-globalDL2.startSeq > 0)
                {
                    xEndLine = Math.round((annots[i].end-globalDL2.startSeq)/globalDL2.scaleSeqScreen);
                }
                else
                {
                    xEndLine = width;
                }
                if(annots[i].sense == "+") yLine=ySenseLine;
                if(annots[i].sense == "-") yLine=yAntisenseLine;

                svg.append("line")
                    .attr("class", globalDL3.classSVG+" line "+annots[i].type)
                    .attr("x1", xStartLine)
                    .attr("x2", xEndLine)
                    .attr("y1", yLine)
                    .attr("y2", yLine);

            }
        }





    }
}



function dataLine2(point, sizePattern, numNucleotides)
{
    // Get information of dataLine2
    var nameSVG=globalDL2.nameSVG;
    var classSVG=globalDL2.classSVG;
    var width=globalDL2.graphWidth;
    var height=globalDL2.graphHeight;
    var x0=globalDL2.paddingX;
    var y0=globalDL2.paddingY;


    // We round to the nearest hundred from sizePattern
    var focusLine = Math.round(sizePattern*globalSeq.ws/100)*100;

    var startSeq = point-(numNucleotides/2)+(focusLine/2);
    var endSeq   = point+(numNucleotides/2)+(focusLine/2);

    // Get the part of sequence that we need to draw
    var result = Server.getPartSeq(startSeq,endSeq);


    var dataLine = drawDataLine(nameSVG, classSVG, width, height, x0, y0,
                                result.partSeq, 1, startSeq, 0, numNucleotides, // scaleSeqServ is 1 because is a part of sequence
                                point, sizePattern);

    // Save information of dataLine2
    globalDL2.width = dataLine.width;
    globalDL2.height = dataLine.height;
    globalDL2.scaleSeqScreen = dataLine.scaleSeqScreen;
    globalDL2.scaleServScreen = dataLine.scaleServScreen;
    globalDL2.data = dataLine.data;
    globalDL2.x = dataLine.x;
    globalDL2.y = dataLine.y;
    globalDL2.svg = dataLine.svg;
    globalDL2.startSeq = startSeq;
    globalDL2.endSeq = endSeq;
}


function drawDataLine(nameSVG, classSVG, width, height, x0, y0, seqServ, scaleSeqServ, initialPoint, startSeq, endSeq, point, sizePattern)
{
    // Default values
    point       || ( point = 0 );
    sizePattern || ( sizePattern = 0 );


    // First, we delete the image, if this exist
    var imageSVG = $("#"+nameSVG);
    if ( imageSVG.length) { imageSVG.empty(); }

    // Sizes...
    var margin = {top: y0, right: x0, bottom: y0, left: x0};
    width = width - margin.left - margin.right;
    height = height - margin.top - margin.bottom;


    // Get info about sequence
    var mean=globalSeq.mean;
    var stdev=globalSeq.stdev;


    // Define the scales with that we will work
    var startData = Math.floor(startSeq/scaleSeqServ);
    var endData = Math.ceil(endSeq/scaleSeqServ);
    var sizeData = endData-startData;
    var sizeSeq = endSeq-startSeq;


    var scaleServScreen = sizeData/width;
    var scaleSeqScreen=1;
    if(sizeSeq>width)   // width/bps compression
    {
        scaleSeqScreen=sizeSeq/width; // i.e. nucleotides per pixel
    }
    if(DEBUG_GBV) console.log("\ndataLine(): startSeq: "+startSeq+" - endSeq: "+endSeq+" - sizeSeq: "+sizeSeq+" - graphWidth: "+width);
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
        .range([0, width]);
    var y = d3.scale.linear()
        .domain([ymin, ymax])
        .range([height, 0]);

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
    svg.append("path")
        .datum(data)
        .attr("class", classSVG+" line")
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
        .attr("x", width-50)
        .attr("y", 15);

    // Save information of dataLine1
    return {
        width : width,
        height : height,
        scaleSeqScreen : scaleSeqScreen,
        scaleServScreen : scaleServScreen,
        data : data,
        x : x,
        y : y,
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
