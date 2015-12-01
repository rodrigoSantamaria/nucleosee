


var globalSeq =
{
    seq : null,
    mean : null,
    stdev : null,
    ws : null,
    scaleSeqServ : null
};

var globalDL1 =
{
    graphHeight : 200,
    graphWidth : screen.width,
    paddingX : 50,
    paddingY : 50,

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
    graphHeight : 200,
    graphWidth : screen.width,
    paddingX : 50,
    paddingY : 50,

    width : null,
    height : null,
    scaleSeqScreen : null,
    scaleServScreen : null,
    data : null,
    x : null,
    y : null,
    svg : null
};


function dataLine1(seq, maxSize, start, end, mean, stdev, ws)
{
    // Get and put global variables...
    globalSeq.seq=seq;
    globalSeq.mean=mean;
    globalSeq.stdev=stdev;
    globalSeq.ws=ws;
    var width=globalDL1.graphWidth;
    var height=globalDL1.graphHeight;
    var x0=globalDL1.paddingX;
    var y0=globalDL1.paddingY;


    // Sizes...
    var margin = {top: y0, right: x0, bottom: y0, left: x0};
    width = width - margin.left - margin.right;
    height = height - margin.top - margin.bottom;
    var sizeSeq = end-start;

    // Define the scales with that we will work
    globalSeq.scaleSeqServ = Math.floor(sizeSeq/maxSize);
    var scaleServScreen = Math.ceil(seq.length/width);
    var scaleSeqScreen=1;
    if(sizeSeq>width)   // width/bps compression
    {
        scaleSeqScreen=Math.floor(sizeSeq/width); // i.e. nucleotides per pixel
    }
    if(DEBUG) console.log("dataLine(): nucleotides/pixel: "+scaleSeqScreen+" - start: "+start+" - end: "+end+" - size: "+sizeSeq+" - graphWidth: "+width);

    // Calculate 'y minimum' and 'y maximum'
    var ymin=0;
    var ymax=mean+3*stdev;
    if(mean-3*stdev>0)
        ymin=mean-3*stdev;

    console.log("scaleSeqScreen: "+scaleSeqScreen+" - (sizeSeq/width) -> ("+sizeSeq+"/"+width+"): "+(sizeSeq/width));
    console.log("scaleSeqServ: "+globalSeq.scaleSeqServ);
    console.log("scaleServScreen: "+scaleServScreen+" - (seq.length/width) -> ("+seq.length+"/"+width+"): "+(seq.length/width));

    // Create the array with a value by pixel
    var data=[];
    for(var i=start,k=0;i<seq.length;i=i+scaleServScreen,k++)
    {
        var average=0;
        for(var j=i;j<i+scaleServScreen;j++)
            if(seq[j]>=0)
                average+=seq[j];
        average/=scaleServScreen;

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
        .tickFormat(function(d) {return roundTickFormat(d,scaleSeqScreen,start,end)})
        .tickValues(getTicks(sizeSeq, scaleSeqScreen, 10));
    var yAxis = d3.svg.axis().scale(y).orient("left");

    // Function to draw the line
    var line = d3.svg.line()
        .x(function(d) { return x(d.pos); })
        .y(function(d) { return y(d.value); });


    // Image SVG: image
    var svg = d3.select("#lineSeq")
        .append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
        .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // Image SVG: axis x
    svg.append("g")
            .attr("class", "dl1 x axis")
        .call(xAxis)
            .attr("transform", "translate(0," + height + ")");

    // Image SVG: axis y
    svg.append("g")
            .attr("class", "dl1 y axis")
        .call(yAxis)
        .append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", 6)
            .attr("dy", ".71em")
            .style("text-anchor", "end");

    // The image SVG: line
    svg.append("path")
        .datum(data)
            .attr("class", "dl1 line")
            .attr("d", line);

    // The image SVG: scale text
    svg.append("text")
            .text("1 : "+scaleSeqScreen)
            .attr("class", "dl1 scale")
            .attr("x", width-50)
            .attr("y", 15);

    // Save information of dataLine1
    globalDL1.width = width;
    globalDL1.height = height;
    globalDL1.scaleSeqScreen = scaleSeqScreen;
    globalDL1.scaleServScreen = scaleServScreen;
    globalDL1.data = data;
    globalDL1.x = x;
    globalDL1.y = y;
    globalDL1.svg = svg;
}


function dataLine1_drawPoints(points)
{

    var seqPoints=[];
    for(var i=0; i<points.length;i++)
    {
        seqPoints.push(points[i] * globalSeq.ws);
    }

    console.log("drawPoints");
    console.log(points);
    console.log(seqPoints);

    var dataPoints=[];
    for(i=0; i<seqPoints.length;i++)
    {
        var dataPoint = Math.floor(seqPoints[i]/globalDL1.scaleSeqScreen);
        dataPoints.push({pos: (globalDL1.data)[dataPoint].pos, value: (globalDL1.data)[dataPoint].value});
    }


    // Mouseover tip
    var tip = d3.tip()
        .attr('class', 'dl1 point-tip')
        .offset([120, 40])
        .html(function(d,i)
        {
            dataLine2(seqPoints[i]);
            return "<strong>" + seqPoints[i] + " position</strong><br>" + Math.round(d.value*100)/100 + " value" + "<br>";
        });

    // Calls tip
    globalDL1.svg.call(tip);

    // Remove all points
    globalDL1.svg.selectAll(".points")
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




function dataLine2(point)
{
    // First, we delete the image, if this exist
    var lineSeq2 = $("#lineSeq2");
    if ( lineSeq2.length) { lineSeq2.empty(); }

    // Get global variables...
    var width=globalDL2.graphWidth;
    var height=globalDL2.graphHeight;
    var x0=globalDL2.paddingX;
    var y0=globalDL2.paddingY;
    var seq = globalSeq.seq;
    var mean = globalSeq.mean;
    var stdev = globalSeq.stdev;


    // Sizes...
    var margin = {top: y0, right: x0, bottom: y0, left: x0};
    width = width - margin.left - margin.right;
    height = height - margin.top - margin.bottom;
    var sizeSeq = width*globalSeq.scaleSeqServ;

    // Define the scales with that we will work
    var scaleServScreen = 1;
    var scaleSeqScreen = 1;
    if(sizeSeq>width)   // width/bps compression
    {
        scaleSeqScreen=Math.floor(sizeSeq/width); // i.e. nucleotides per pixel // TODO: quitar el math.floor
    }
    if(DEBUG) console.log("dataLine2(): nucleotides/pixel: "+scaleSeqScreen+" - size: "+sizeSeq+" - graphWidth: "+width);


    // Calculate 'y minimum' and 'y maximum'
    var ymin=0;
    var ymax=mean+3*stdev;
    if(mean-3*stdev>0)
        ymin=mean-3*stdev;


    // Create the array with a value by pixel
    var data=[];
    //var startData = Math.floor(point/24)-(width/2);
    //var endData = Math.floor(point/24)+(width/2)-1;
    var numNucleotides = 2000;
    var startData = Math.floor(point/24)-Math.floor(numNucleotides/2/24);
    var endData = Math.floor(point/24)+Math.floor(numNucleotides/2/24)-1;
    for(var i=startData,k=0;i<endData;i=i+scaleServScreen,k++)
    {
        var average=0;
        for(var j=i;j<i+scaleServScreen;j++)
            if(seq[j]>=0)
                average+=seq[j];
        average/=scaleServScreen;

        if (average<ymin) average=ymin;
        if (average>ymax) average=ymax;

        data.push({pos: k, value: average})
    }
    console.log(data);


    // Scaling of the axes
    var x = d3.scale.linear()
        .domain(d3.extent(data, function (d) { return d.pos; }))  // xmin, xmax
        .range([0, width]);
    var y = d3.scale.linear()
        .domain([ymin, ymax])
        .range([height, 0]);

    // Axis labels
    var xAxis = d3.svg.axis().scale(x).orient("bottom")
        .tickFormat(function(d) {return roundTickFormat(d,scaleSeqScreen,startData*globalSeq.scaleSeqServ,endData*globalSeq.scaleSeqServ)})
        .tickValues(getTicks(sizeSeq, scaleSeqScreen, 10));
    var yAxis = d3.svg.axis().scale(y).orient("left");

    // Function to draw the line
    var line = d3.svg.line()
        .x(function(d) { return x(d.pos); })
        .y(function(d) { return y(d.value); });


    // Image SVG: image
    var svg = d3.select("#lineSeq2")
        .append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
        .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // Image SVG: axis x
    svg.append("g")
            .attr("class", "dl2 x axis")
        .call(xAxis)
            .attr("transform", "translate(0," + height + ")");

    // Image SVG: axis y
    svg.append("g")
            .attr("class", "dl2 y axis")
        .call(yAxis)
        .append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", 6)
            .attr("dy", ".71em")
            .style("text-anchor", "end");

    // The image SVG: line
    svg.append("path")
        .datum(data)
            .attr("class", "dl2 line")
            .attr("d", line);

    // The image SVG: scale text
    svg.append("text")
        .text("1 : "+scaleSeqScreen)
            .attr("class", "dl2 scale")
            .attr("x", width-50)
            .attr("y", 15);

    // Save information of dataLine1
    globalDL2.width = width;
    globalDL2.height = height;
    globalDL2.scaleSeqScreen = scaleSeqScreen;
    globalDL2.scaleServScreen = scaleServScreen;
    globalDL2.data = data;
    globalDL2.x = x;
    globalDL2.y = y;
    globalDL2.svg = svg;
}







function getTicks(sizeSeq, window)
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
        ticks.push((factorLabel+i*factorLabel)/window);
    }

    return ticks;
}



function roundTickFormat(d, window, start, end)
{
    var tickLabel = start+d*window;

    if (tickLabel != "0")
    {
        if(tickLabel/1000000 >= 1 && (end-start)>10000000)
            tickLabel=Math.round(tickLabel/1000000)+"M";
        else if(tickLabel/1000>=1 && (end-start)>10000)
            tickLabel=Math.round(tickLabel/1000)+"K";
    }

    return tickLabel;
}









