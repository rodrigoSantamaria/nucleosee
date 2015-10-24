/**
 * Created by jonatan on 15/10/15.
 */




function dataLine(seq, start, end, mean, stdev, height, width, x0, y0)
{

    //Compute zoom ratio
    var window=1;
    var xdist=width/(end-start);

    if((end-start)>width)       // width/bps compression
    {
        window=Math.ceil((end-start)/width); // i.e. nucleotides per pixel
        xdist=1;
    }

    // Sizes...
    var margin = {top: y0, right: x0, bottom: y0, left: x0},
        width = width - margin.left - margin.right,
        height = height;

    console.log("dataLine(): nucleotides/pixel: "+window+" - start: "+start+" - end: "+end+" - size: "+(end-start)+" - graphWidth: "+width);


    // Scaling of the axes
    var x = d3.scale.linear()
        .range([0, width]);

    var y = d3.scale.linear()
        .range([height, 0]);

    // Axis labels
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    // Function to draw the line
    var line = d3.svg.line()
        .x(function(d) { return x(d.pos); })
        .y(function(d) { return y(d.value); });

    // The image SVG
    var svg = d3.select("#dna")
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


    ymax=mean+3*(stdev/2);
    if(mean-3*(stdev/2)>0)
        ymin=mean-3*(stdev/2);
    else
        ymin=0;


    data= [];
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

    console.log("ymin:"+ymin+" ymax:"+ymax);
    x.domain(d3.extent(data, function (d) { return d.pos; }));
    y.domain([ymin, ymax]);



    //Mouseover tip
    var tip = d3.tip()
        .attr('class', 'd3-tip')
        .offset([120, 40])
        .html(function(d) {
            dataLine2(seq, (d.pos*window)-(graphWidth/2), (d.pos*window)+(graphWidth/2), mean, stdev, height, width, x0, y0);
            return "<strong>" + d.pos + " position</strong><br>" + d.value + " value" + "<br>";
        });

    svg.call(tip);


    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")");

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis)
        .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end")

    svg.append("path")
        .datum(data)
        .attr("class", "line")
        .attr("d", line);


    data2=[];
    data2.push({pos: data[900].pos, value: data[900].value})
    data2.push({pos: data[600].pos, value: data[600].value})
    data2.push({pos: data[300].pos, value: data[300].value})


    svg.selectAll(".data2")
        .data(data2)
        .enter().append("circle")
        .attr('class', 'datapoint')
        .attr('cx', function(d) { return x(d.pos); })
        .attr('cy', function(d) { return y(d.value); })
        .attr('r', 6)
        .attr('fill', 'white')
        .attr('stroke', 'steelblue')
        .attr('stroke-width', '3')
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);

}


function dataLine2(seq, start, end, mean, stdev, height, width, x0, y0)
{

    if ( $("#dna2").length)  { $("#dna2").empty(); }


    //Compute zoom ratio
    var window=1;
    var xdist=width/(end-start);

    if((end-start)>width)       // width/bps compression
    {
        window=Math.ceil((end-start)/width); // i.e. nucleotides per pixel
        xdist=1;
    }
    console.log("dataLine(): nucleotides/pixel: "+window+" - start: "+start+" - end: "+end+" - size: "+(end-start)+" - graphWidth: "+width);



    // Sizes...
    var margin = {top: y0, right: x0, bottom: y0, left: x0},
        width = width - margin.left - margin.right,
        height = height;

    // Scaling of the axes
    var x = d3.scale.linear()
        .range([0, width]);

    var y = d3.scale.linear()
        .range([height, 0]);

    // Axis labels
    var xAxis = d3.svg.axis()
        .scale(x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

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


    ymax=mean+3*(stdev/2);
    if(mean-3*(stdev/2)>0)
        ymin=mean-3*(stdev/2);
    else
        ymin=0;


    data= [];
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

    x.domain(d3.extent(data, function (d) { return d.pos; }));
    y.domain([ymin, ymax]);



    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")");

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis)
        .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end")

    svg.append("path")
        .datum(data)
        .attr("class", "line")
        .attr("d", line);

}