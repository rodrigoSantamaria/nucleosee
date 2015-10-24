// Creates canvas 320 × 200 at 10, 50
screenHeight=500;
screenWidth=1200;
startX=20;
startY=50;
var paper = new Raphael(0, startY, screenWidth, screenHeight);
var hist;


/**
 * Draws a histogram given a sequence w
 * @param w sequence of real numbers
 * @param max max value of w (or smaller if wanted). This is done to avoid recomputation (although possible unnecessary)
 * @param min min value of w (or higher if wanted). This is done to avoid recomputation (although possible unnecessary)
 * @param numBins numer of bins to divide the histogram into
 * @param x0 position in the screen of top-left corner of the histogram
 * @param y0
 * @param histHeight height of the histogram
 * @param color
 * @param opacity
 */
function histogram(w, max, min, numBins, x0, y0, histHeight, color, opacity)
    {
    this.freqs=distribution(w,max, min,numBins);
    this.maxf=maximum(this.freqs);
    var p;
    for(i in this.freqs)
        p+="M"+(x0+parseFloat(i))+","+y0+"l"+0+","+(0-(this.freqs[i]/this.maxf*histHeight));
    this.visual=paper.path(p);
    this.visual.attr({"fill": color});
    this.visual.attr({"stroke": color});
    this.visual.attr({"opacity": opacity});
    }

function sinusoid(w, x0, y0, height)
{
    var xi=x0;
    var yi=y0+Math.sin(x0)
    var p="M"+xi+","+yi;
    for(var i=0;i< w.length;i++)
    {
        var yn=Math.sin(x0+i);
        p+="L"+(x0+i)+","+(y0+height*yn);
        yi=yn;
    }
    this.visual=paper.path(p);
}

/**
 * Draws a segment with the given info

 * @param line  Visual element over which to get a segment
 * @param segmentStart starting position of line from which to draw
 * @param segmentEnd
 * @param disc
 * @param yDesp
 * @param tails
 * @param color
 * @param opacity
 *//*
 function segment(line, segmentStart, segmentEnd, disc, yDesp, tails, color, opacity)
 {
     var sp=line.visual.getPath().slice(Math.round(segmentStart*screenWidth/disc.length), Math.round(segmentEnd*screenWidth/disc.length));
     if(sp.length>0)
         {
         sp[0][0]="M";
         for(i=0;i<sp.length;i++)
             {
                 if(sp[i][2]<100)           //super horrid hack to make it work on several calls... why is it persistent???
                     sp[i][2]+=yDesp;
             }
         this.visual=paper.path(sp);
         this.visual.attr({"stroke":color, "opacity":opacity});
         }
     else
        {
        this.visual=paper.path("M0,0L1,1");
        this.visual.attr({"stroke":"#000000", "opacity":0});
        }
 }*/

/**
 * Creates a visual element with the sequence especified at na. Notice this is a detailed line, rarely drawn
 * completely (as it is both useless, costly and impossible to draw on a normal screen).
 * However, once the visual element is build, we can use method SEGMENT to draw small pieces on demand
 * @param na    sequence with all the values (typically a chromosome or other large sequence)
 * @param start
 * @param end
 * @param max
 * @param min
 * @param height height and width of the graph in which to draw the line
 * @param width
 * @param x0     x0,y0 of the top-left square of the graph in which to draw the line
 * @param y0
 * @param linewidth width for the visual line
 * @param color
 * @param opacity
 * @param animationTime
 */



function dataLine(na, start, end, max, min, height, width, x0, y0, linewidth, color, opacity, animationTime, drawScale)
{
    //Compute zoom ratio
    var window=1;
    var xdist=width/(end-start);
    if((end-start)>width)       // width/bps compression
        {
        window=Math.ceil((end-start)/width); //i.e. nucleotides per pixel
        xdist=1;
        }
    console.log("   --> dataLine: window: "+window+"\tstart: "+start+"\tend: "+end+"\tsize: "+(end-start)+"\tgraphWidth: "+width);
    //console.log("x0: "+x0);
    //INITIAL DATA LINE
    var p="M 10 50";
    for(i=0;i<100;i++)
        {
        p+=" l "+screenWidth/100+" "+y0;
        }
    this.visual=paper.path(p);
    if(drawScale) {
        this.visual.scale = paper.text(x0 + width - 20, y0 + 10, "1:" + window);
        this.visual.scale.attr({'stroke': color, 'opacity': opacity});
        }

    //DRAW DATA LINE -----------------------------------------------------
    var yant=y0+height-height*(na[start]-min)/(max-min);
    //console.log("max y min: "+max+", "+min);
    var points=new Array();
    var p="M "+x0+","+yant;
    points.push(p);

    for(i=start;i<end;i+=window)
        {
        var average=0;
        for(j=i;j<i+window;j++)
            average+=na[j];
        average/=window;

        var y=y0+height-height*(average-min)/(max-min);
        var s=" l "+xdist+","+(y-yant);
        p+=s
        points.push(s);
        yant=y;
        }
    //if(window==1)
    //    console.log("path: "+p);
    this.visual.animate({path: p}, animationTime, '>');
    this.data=na;
    this.data.rawPath=p;
    this.visual.attr({'stroke-width':linewidth, 'stroke':color, 'opacity':opacity});
    }

/**
 * Draws a horizontal line with tickmarks
 * @param fullLength length of the sequence that the ticks respond to (full or original, normalized sequence)
 * @param x0 starting position of the graphic, top left (line will be drawn at its bottom)
 * @param y0
 * @param height height of the graphic (so we can compute where's its bottom)
 * @param width width of the graphic (so we can compute where's its bottom)
 * @param startTick first tick value to draw (usually 0)
 * @param numTicks number of ticks to draw after the first one (approximate, the algorithm may adjust a little)
 * @param color
 * @param opacity
 * @param mean  #TODO: for vertical ticks (by now only horizontal)
 * @param max
 * @param min
 */
function ticks(fullLength, x0, y0, height, width, startTick, numTicks, color, opacity, mean, max, min)
    {
    var numZeroes=-0
    while(fullLength/Math.pow(10,numZeroes)>numTicks)
        {
        numZeroes += 1
        }
    if(fullLength/Math.pow(10,numZeroes)<4)
        numZeroes-=1
    var factorLabel=Math.pow(10, numZeroes);
    this.labels=[];
    this.marks=[];
    numTicks=Math.round(fullLength/Math.pow(10,numZeroes))
    console.log("   --> ticks: fullLength: "+fullLength+"\tnumZeroes: "+numZeroes+"\tnumTicks: "+numTicks);
    for(var i=0;i<numTicks;i++)
        {
        var y=y0+10+height;
        var x=x0+i*width*factorLabel/fullLength;

            var tickLabel=""
        if(i==0)    tickLabel=startTick
        else        tickLabel=startTick-startTick%factorLabel+i*factorLabel;
        if(tickLabel/1000000 >= 1)
            tickLabel=tickLabel/1000000+"M"
        else if(tickLabel/1000>=1)
            tickLabel=tickLabel/1000+"K"
        this.labels[i]=paper.text(x,y,tickLabel);
        this.labels[i].attr({'stroke':color, 'opacity':opacity})

        this.marks[i]=paper.path("M"+x+","+(y0+height-3)+"L"+x+","+(y0+height+3));
        this.marks[i].attr({'stroke':color, 'opacity':opacity})
        }

    //final tick (provides the full length, rounded
    var x=(x0+width)
    tickLabel=startTick+fullLength;
    if(tickLabel/1000000 >= 1)
        tickLabel=Math.round(tickLabel/1000000)+"M"
    else if(tickLabel/1000>=1)
        tickLabel=Math.round(tickLabel/1000)+"K"
    this.labels[numTicks]=paper.text(x,y,tickLabel);

    this.labels[numTicks].attr({'stroke':color, 'opacity':opacity})
    this.marks[numTicks]=paper.path("M"+x+","+(y0-3)+"L"+x+","+(y0+3));
    this.marks[numTicks].attr({'stroke':color, 'opacity':opacity})

    //GUIDE LINES
    var p="M "+x0+" "+(y0+height-height*(mean-min)/(max-min));
    for(i=0;i<100;i++)
        p+=" l "+width/100+" "+0;
    //console.log("mean line at "+(height-height*(mean-min)/(max-min)));
    this.meanLine=paper.path(p);
    this.meanLine.attr(
        {
            stroke: '#111111',
            'stroke-width': 0.1,
            'stroke-dasharray' : "·"
        }
    );
    }

/**
 * Computes portions of the sequence that do not vary on the discrete count
 * For example, if we have a-a-a-a-a-b-b-c-a-a-a-a we will have four segments
 * @param x0
 * @param graphHeight
 * @param graphWidth
 * @param disc
 * @param numBins
 * @param ppp numerical points represented per pixel
 * @param width stroke characteristics
 * @param color
 * @param opacity
 */
function segmentTrack(x0, graphHeight,graphWidth,disc, numBins, ppp, width, color, opacity)
    {
    var interval=graphHeight/5;
    var j=0;
    this.segments=[]
    this.width=0

    for(i=0;i<disc.length-1;i++)
        {
        var pos=disc[i].charCodeAt()-'a'.charCodeAt();
        var y=graphHeight-graphHeight*(pos/numBins);
        var nextPos=disc[i+1].charCodeAt()-'a'.charCodeAt();
        var w=1;
        var startSegment=i;
        while(nextPos==pos && i+1<disc.length)
            {
            i++;
            w++
            nextPos=disc[i].charCodeAt()-'a'.charCodeAt();
            }
        if(w>1)
            {
            i--;
            w--;
            }

        p="M "+(x0+startSegment*(graphWidth/disc.length))+","+y+" l "+w*graphWidth/disc.length+", "+0;
        this.segments[j]=paper.path(p);
        this.segments[j].attr({'stroke-width':width, 'stroke':color, 'opacity':opacity});

        //visual start and size
        this.segments[j].width=w*graphWidth/disc.length; //visual width in pixels
        this.segments[j].x0=startSegment*(graphWidth/disc.length);//visual start (respect to segmentTrack)

        //numeric start and size
        this.segments[j].start=startSegment*ppp;        //numeric start respect to na
        this.segments[j].end=startSegment*ppp+w*ppp;    //end

        this.segments[j].letter=disc[startSegment];    //discretization letter for the
        this.segments[j].status="normal"               //normal, hovered, clicked
        /*this.segments[j].onclick=function()
            {
            if(this.status=="normal")
                {
                this.attr({'stroke': '#00FF00'});
                this.status="clicked";
                }
            else
                this.attr({'stroke':color});
            }
        this.segments[j].onmouseover=function() {
            this.attr({'stroke': '#FF0000'});
            }
        this.segments[j].onmouseout=function() {
                this.attr({'stroke': color});
            }*///This not working, not sure why

        this.width+=this.segments[j].width;             //full width of the track
        j++;
        }
    }



//************************************************************************
//************************************************************************
//************************************************************************
// ---------------------------- Main code --------------------------------
//************************************************************************
//************************************************************************
//************************************************************************

if (window.File && window.FileReader && window.FileList && window.Blob) //check if there's support for file reading
    {
    console.log("Start all...");
    this.handleEvent=function(evt)
        {
        switch(event.type)
            {
            case 'change':
                console.log("----- UPLOAD FILE -----");

                var files = evt.target.files; // FileList object
                f=files[0]  //By now, just one file

                var output = [];
                graphHeight=100;
                graphWidth=screenWidth-startX*2
                console.log("Name of file: "+f.name);

                //Remove previous visual elements if any (reset)
                line=rawSegment=rawSegmentContext=strack=segmentTicks=lineTicks=undefined
                histSegment=lane=line=undefined


                //0) READ FILE
                //var fr=new FileReader(f); //on the client side
                //var freqs=fr.readAsText(f);
                //var lines=fr.result.split("\n").slice(2);
                //console.log("File with "+lines.length+" lines");
                startTime=new Date()
                var serverFilePath=sendFile(f);    //passing it to the server side (best solution for >1MB files)
                console.log("Path in remote: "+serverFilePath);
                console.log("Time spent sending: "+ (new Date()-startTime)+"ms");


                var x=10;
                var ws=150; //window size: discrete to real ratio
                var nb=5; //num bins
                var maxSize=100000 //maximum number of normalized data to store

                //PREPROCESSING; NORMALIZE, STATISTICAL DESCRIPTORS, DISCRETIZE ----------------------------------------
                startTime=new Date()
                var desc=preprocess(serverFilePath, ws, nb, maxSize);    //normalize+stats+discretize+suffix array build (next)

                max=desc.max;
                min=desc.min;
                fullLength=desc.fullLength
                var na=desc.seq;    //just a sampling of about 100K of the original full length sequence
                var disc=desc.dseq; //whole seq compression


                var p="";
                console.log("Time spent preprocessing: "+ (new Date()-startTime)+"ms");

                console.log("Length of na:"+na.length+" (full length="+desc.fullLength+")");


                //send + process a full Sc wig is about 70s+50s ~2 minutes.
                    //1) Do not resend if already uploaded
                    //2) Either save processing to not re-preprocess, and/or think about binary wigs

                var i=1
                //DRAW BACKGROUND LANE
                console.log("----- DRAWING BACKGROUND LANE -----");
                lane=paper.rect(startX, 0+graphHeight*(i-1), graphWidth, graphHeight); //for some reason, i is already increased!
                lane.attr({'fill':'#FFFFFF', 'stroke':'#AAAAAA', 'fill-opacity':'0'});

                //BACGROUND LANE FOR 2nd TRACK
                lane2=paper.rect(startX, 20+graphHeight*(i), graphWidth, graphHeight); //for some reason, i is already increased!
                lane2.attr({'fill':'#FFFFFF', 'stroke':'#AAAAAA', 'fill-opacity':'0'});

                //DATA LINE (preprocessed data)
                line=new dataLine(na,0, na.length,max, min, graphHeight, graphWidth, startX,0, '1', '#000000', '0.2', 1000, true);

                //TICK MARKS
                lineTicks=new ticks(desc.fullLength, startX, 0, graphHeight, graphWidth, 0, 10, '#000000', '0.2', desc.mean, desc.max, desc.min);

                //DRAW DISCRETE SEGMENTS (discretization by letter)
                strack=new segmentTrack(startX, graphHeight,graphWidth,disc,nb, ws, '5', '#0000FF', '0.2');

                //DATA LINE INTERACTION ----------------------------------
                lineSegment=histSegment=rawSegment=undefined;
                var hseg=-1; //hovered segment



                var hoverIn=function(e)
                    {
                    if(histSegment!=undefined)  hoverOut(e);

                    //discretized line segment
                    for(var j=0;j<strack.segments.length;j++)
                        {
                        if(startX+strack.segments[j].x0 < e.offsetX && (startX+strack.segments[j].x0+strack.segments[j].width) > e.offsetX)
                            {
                            hseg=j;
                            strack.segments[j].attr({'stroke':'#FF0000'});
                            break
                            }
                        }
                    var segmentStart=strack.segments[hseg].x0; //visual start
                    var segmentEnd=strack.segments[hseg].x0+strack.segments[hseg].width;//visual end

                    var segmentWidth=strack.segments[hseg].width;//visual width

                    //draw context in case it is smaller
                    var segmentStartContext= strack.segments[hseg].start;
                    var segmentEndContext= strack.segments[hseg].end;

                    var rsWidth=strack.segments[hseg].end-strack.segments[hseg].start;
                    var segmentWindowSize=rsWidth/graphWidth;

                     //To avoid under-information on small segments
                     if(segmentWindowSize<1)
                         {
                         var freeSpace=graphWidth*(1-segmentWindowSize);
                         segmentStartContext=Math.max(0,strack.segments[hseg].start-0.5*freeSpace);
                         segmentEndContext=Math.min(strack.segments[hseg].end+0.5*freeSpace, na.length);
                         console.log("FreeSpace: "+freeSpace+"\tgW: "+graphWidth+"\tsW: "+segmentWidth);
                         }

                    console.log("expanded segment");
                    //context data line (black - dark red)
                    rawSegmentContext=new dataLine(na,Math.round(segmentStartContext), Math.round(segmentEndContext),max, min, graphHeight, graphWidth, startX,graphHeight+20, '1', '#330000', '1', 0, true);
                    //selected segment data line (red)
                    console.log("real segment");
                    var esWidth=segmentEndContext-segmentStartContext;

                    //TODO: segments in 1:2+ do not show the rigth ticks (remaining fragment not drawn...)
                    //TODO: last segments are incorrectly drawn (displaced)
                    rawSegment=new dataLine(na,strack.segments[hseg].start, strack.segments[hseg].end,max, min, graphHeight, graphWidth*(rsWidth/esWidth), startX+(strack.segments[hseg].start-segmentStartContext),graphHeight+20, '2', '#FF0000', '1', 0, false);
                    segmentTicks=new ticks(Math.round(segmentEndContext)-Math.round(segmentStartContext), startX, graphHeight+20, graphHeight, graphWidth, segmentStartContext,10, '#880000', '0.2', desc.mean, desc.max, desc.min);
                    var nas=na.slice(strack.segments[hseg].start, strack.segments[hseg].end)
                    histSegment=new histogram(nas, max, min, 50, 15, 40, 30, "#FF0000", 0.5);
                    }

                var hoverOut=function(e)    //clear the histogram/line
                    {
                    //if(rawSegment.visual.scale!=null)    rawSegment.visual.scale=undefined;
                    //if(rawSegmentContext.visual.scale!=null)        rawSegmentContext.visual.scale=undefined;
                    rawSegment.visual.remove();
                    rawSegmentContext.visual.remove();

                        if(hseg>=0)    strack.segments[hseg].attr({'stroke':'#0000FF'});
                    for(var i=0;i<segmentTicks.marks.length;i++)
                        {
                        segmentTicks.marks[i].remove();
                        segmentTicks.labels[i].remove();
                        }
                    segmentTicks.meanLine.remove();
                    //segmentTicks.zoomLabel.remove();

                    histSegment.visual.remove();
                    }

                var clickMethod=function()
                    {
                        console.log("here click");
                    }

                lane.hover(hoverIn, hoverOut, paper, paper);
                lane.click(clickMethod, paper);
                lane.mousemove(hoverIn);
                lane.mouseover(hoverIn);
                lane.mouseout(hoverOut);
                lane.toFront();

                //DRAW DATA DISTRIBUTION
                hist= new histogram(na, maximum(na), minimum(na), 50, 15, 40, 30, "#999999", 1);


                //DRAW sinusoidal
                //sinus=new sinusoid(na, 10, 200, 10);
            }
        }
    document.getElementById('files').addEventListener('change', this, false); //TODO: changing to a post method to upload files to server
    }
else
    {
    alert('The File APIs are not fully supported by your browser.');
    }


