


var globalSeq =
{
    seq : null,
    mean : null,
    stdev : null,
    ws : null
};

var globalDL1 =
{
    graphHeight : 200,
    graphWidth : screen.width,
    startX : 50,
    startY : 50,
    svg : null,
    data : null,
    window : null,
    x : null,
    y : null
};




function dataLine1(seq, start, end, mean, stdev, ws)
{
    globalSeq.seq=seq;
    globalSeq.mean=mean;
    globalSeq.stdev=stdev;
    globalSeq.ws=ws;
    var width=globalDL1.graphWidth;
    var height=globalDL1.graphHeight;
    var x0=globalDL1.startX;
    var y0=globalDL1.startY;


    // Sizes...
    var margin = {top: y0, right: x0, bottom: y0, left: x0};
    width = width - margin.left - margin.right;
    height = height - margin.top - margin.bottom;
    var sizeSeq = end-start;

    // Compute zoom ratio
    var window=1;
    var xdist=width/sizeSeq;

    if(sizeSeq>width)   // width/bps compression
    {
        window=Math.ceil(sizeSeq/width); // i.e. nucleotides per pixel
        xdist=1;
    }
    if(DEBUG) console.log("dataLine(): nucleotides/pixel: "+window+" - start: "+start+" - end: "+end+" - size: "+sizeSeq+" - graphWidth: "+width);


    // Calculate 'y minimum' and 'y maximum'
    var ymin=0;
    var ymax=mean+3*stdev;
    if(mean-3*stdev>0)
        ymin=mean-3*stdev;

    // Create the array with a value by pixel
    var data=[];
    for(i=start,k=0;i<end;i+=window,k++)
    {
        var average=0;
        for(j=i;j<i+window;j++)
            if(seq[j]>0)
                average+=seq[j];
        average/=window;

        if (average<ymin) average=ymin;
        if (average>ymax) average=ymax;

        data.push({pos: k, value: average})
    }


    // Scaling of the axes
    var x = d3.scale.linear().range([0, width]).domain(d3.extent(data, function (d) { return d.pos; }));
    var y = d3.scale.linear().range([height, 0]).domain([ymin, ymax]);

    // Axis labels
    var xAxis = d3.svg.axis().scale(x).tickFormat(function(d) {return roundTickFormat(d,window,start,end)}).orient("bottom").tickValues(getTicks(sizeSeq, window, 10));
    var yAxis = d3.svg.axis().scale(y).orient("left");


    // Function to draw the line
    var line = d3.svg.line()
        .x(function(d) { return x(d.pos); })
        .y(function(d) { return y(d.value); });


    // The image SVG: image
    var svg = d3.select("#lineSeq")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


    // The image SVG: axis x
    svg.append("g")
        .attr("class", "x axis")
        .call(xAxis)
        .attr("transform", "translate(0," + height + ")");

    // The image SVG: axis y
    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis)
        .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end")

    // The image SVG: line
    svg.append("path")
        .datum(data)
        .attr("class", "line")
        .attr("d", line);


    globalDL1.svg = svg;
    globalDL1.data = data;
    globalDL1.x = x;
    globalDL1.y = y;
    globalDL1.window = window;
}


function dataLine2(DEBUG, seq, start, end, mean, stdev, height, width, x0, y0)
{

    // First, we delete the image, if this exist
    if ( $("#lineSeq2").length)  { $("#lineSeq2").empty(); }

    // Sizes...
    var margin = {top: y0, right: x0, bottom: y0, left: x0};
    width = width - margin.left - margin.right;
    height = height - margin.top - margin.bottom;
    var sizeSeq = end-start;

    // Compute zoom ratio
    var window=1;
    var xdist=width/sizeSeq;

    if(sizeSeq>width)   // width/bps compression
    {
        window=Math.ceil(sizeSeq/width); // i.e. nucleotides per pixel
        xdist=1;
    }
    if(DEBUG) console.log("dataLine2(): nucleotides/pixel: "+window+" - start: "+start+" - end: "+end+" - size: "+sizeSeq+" - graphWidth: "+width);


    // Calculate 'y minimum' and 'y maximum'
    var ymin=0;
    var ymax=mean+3*stdev;
    if(mean-3*stdev>0)
        ymin=mean-3*stdev;

    // Create the array with a value by pixel
    var data=[];
    for(i=start,k=0;i<end;i+=window,k++)
    {
        var average=0;
        for(j=i;j<i+window;j++)
            if(seq[j]>0)
                average+=seq[j];
        average/=window;

        if (average<ymin) average=ymin;
        if (average>ymax) average=ymax;

        data.push({pos: k, value: average})
    }


    // Scaling of the axes
    var x = d3.scale.linear().range([0, width]).domain(d3.extent(data, function (d) { return d.pos; }));
    var y = d3.scale.linear().range([height, 0]).domain([ymin, ymax]);


    // Axis labels
    var xAxis = d3.svg.axis().scale(x).tickFormat(function(d) {return roundTickFormat(d,window,start,end)}).orient("bottom").tickValues(getTicks(sizeSeq, window, 10));
    var yAxis = d3.svg.axis().scale(y).orient("left");


    // Function to draw the line
    var line = d3.svg.line()
        .x(function(d) { return x(d.pos); })
        .y(function(d) { return y(d.value); });


    // The image SVG
    var svg = d3.select("#lineSeq2")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // The image SVG: axis x
    svg.append("g")
        .attr("class", "x axis")
        .call(xAxis)
        .attr("transform", "translate(0," + height + ")");

    // The image SVG: axis y
    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis)
        .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end")

    // The image SVG: line
    svg.append("path")
        .datum(data)
        .attr("class", "line")
        .attr("d", line);
}



function drawPoints(points)
{
    console.log(points);

    var dataPoints=[];
    for(var i=0; i<points.length;i++)
    {
        var point = Math.ceil(points[i] * globalSeq.ws/globalDL1.window);
        dataPoints.push({pos: (globalDL1.data)[point].pos, value: (globalDL1.data)[point].value});
    }


    // Mouseover tip
    var tip = d3.tip()
        .attr('class', 'd3-tip')
        .offset([120, 40])
        .html(function(d)
        {
            var algo = (screen.width/2)-50;
            dataLine2(true, globalSeq.seq, (d.pos*globalDL1.window)-algo, (d.pos*globalDL1.window)+algo, globalSeq.mean, globalSeq.stdev, 200, screen.width, 50, 50);
            return "<strong>" + d.pos*globalDL1.window + " position</strong><br>" + Math.round(d.value*100)/100 + " value" + "<br>";
        });


    globalDL1.svg.call(tip);

    globalDL1.svg.selectAll("circle.datapoint").remove();

    globalDL1.svg.selectAll(".data")
        .data(dataPoints)
        .enter().append("circle")
        .attr('class', 'datapoint')
        .attr('cx', function(d) { return globalDL1.x(d.pos); })
        .attr('cy', function(d) { return globalDL1.y(d.value); })
        .attr('r', 6)
        .attr('fill', 'white')
        .attr('stroke', 'steelblue')
        .attr('stroke-width', '3')
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);
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









