/**
 * Created by jonatan on 15/10/15.
 */




function dataLine(DEBUG, seq, start, end, mean, stdev, height, width, x0, y0)
{
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
    var svg = d3.select("#dna")
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


    // We refund the information necessary to put points
    var result=[];
    result.svg=svg;
    result.data=data;
    result.x=x;
    result.y=y;
    result.window=window;
    return result;
}


function dataLine2(DEBUG, seq, start, end, mean, stdev, height, width, x0, y0)
{

    // First, we delete the image, if this exist
    if ( $("#dna2").length)  { $("#dna2").empty(); }

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
    var svg = d3.select("#dna2")
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



function drawPoints(DEBUG, DL1_svg, data, DL1_x, DL1_y, DL1_window,
                     seq, mean, stdev, height, width, x0, y0){


    // Mouseover tip
    var tip = d3.tip()
        .attr('class', 'd3-tip')
        .offset([120, 40])
        .html(function(d)
        {
            var algo = (width/2)-x0;
            dataLine2(true, seq, (d.pos*DL1_window)-algo, (d.pos*DL1_window)+algo, mean, stdev, height, width, x0, y0);
            return "<strong>" + d.pos*DL1_window + " position</strong><br>" + Math.round(d.value*100)/100 + " value" + "<br>";
        });

    DL1_svg.call(tip);


    DL1_svg.selectAll(".data")
        .data(data)
        .enter().append("circle")
        .attr('class', 'datapoint')
        .attr('cx', function(d) { return DL1_x(d.pos); })
        .attr('cy', function(d) { return DL1_y(d.value); })
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



    //for(var i=0;i<numTicks;i++)------------------------------------------------------------------
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









