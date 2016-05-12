/*
  ┌────────────────────────────────────────────────────────────┐
  │ main.js                                                    │
  ├────────────────────────────────────────────────────────────┤
  │ Description:                                               │
  └────────────────────────────────────────────────────────────┘
*/


var user="jpiriz";
var password="ninguna";

var DEBUG_GBV = false;

/*
 ┌────────────────────────────────────────────────────────────┐
 │                          TODO list                         │
 └────────────────────────────────────────────────────────────┘
  Frontend
  - Permitir sobreimponer sobre DL1 una parrilla explicativa del binarizado, por ejemplo con franjas
  de colores alternos indicando los rangos de a,b,c,d,e, etc. y quizas tambien el tamaño de la ventana
  - Selección de cromosomas
  - Reloj/barra de espera para cargas y búsquedas
  - Dejar el color interno de un punto al hacer mouseout (quitarlo solo al hacer un nuevo mousein)
  - Estilos de la web (sobre todo de la parte de arriba)
  - Pantalla de login

  Backend
  - Intentar optimizar el cálculo de estadísticas y discretización (2.5/5s en 2M y  5/7s en 4M)
  - Intentar optimizar los tiempos de lectura y formateo (en total unos 10s para pombe entero)
        DONE: Los dos puntos anteriores están bastante optimizados. Para pombe entero estamos hablando de 10s para todo
            el proceso ahora, incluyendo lectura, formateo, estadísticas, discretización y bwt.
  - Guardar datos para mayor optimización en cargas posteriores
        DONE: mediante pickle y con un retardo de 2s en la carga inicial (de 12 a 14s) tenemos luego unas cargas de 4s o menos
        TODO: explorar pandas como alternativa a pickle, sobre todo pensando en humano
        TODO: también pensando en humano se puede pensar en integrar el BWT u otros cálculos en el fichero que se guarde
        TODO: se puede también pensar en NO guardar el .wig original, sólo nuestros datos, para evitar duplicaciones de memoria en disco
  - Soporte a la selección de cromosomas
        En el servidor con los cambios necesarios en cliente.
        Básicamente, las interfaces incluyen normalmente un argumento nuevo, 'track', que de momento lo estamos poniendo
        a processedData.chromosomes[0] pero que luego cambiará en función de un combo box.
  - Crear carpetas y fichero de contraseñas + comprobaciones
  - Gestión de distintos organismos (ver cómo diseñarlo)
  - Dar soporte a bigwig
  - Enriquecimiento: estamos teniendo en cuenta unos GIs que se han escogido en un cromosoma pero usando como universo
        el genoma entero. O bien debemos tomar todos los GIs del genoma (quizás sería lo suyo, aunque habría que ver
        cómo responde) o bien tomar como universo sólo los genes de ese cromosoma. De momento lo que hacemos es poner
        un FDR muy estricto, pero esto no es más que un parche.
        Ejemplo: a*5+abcba(0) en ch2 de pombe da 281 ocurrencias de las cuales 39 están en genes. El universo completo
        son 5346 genes con alguna anotación GO --> un FDR de 0.01 no da todavía muchas anotaciones relevantes!

        La idea entonces sería tener preprocesados TODOS los cromosomas? En pombe serían unos 3s por cromosoma, es
        decir un tiempo adicional de 6s, pasando de 14 a 20s de preprocesamiento. En organismos grandes es inasumible
        yo creo. También multiplica por el número de cromosomas el tiempo de la búsqueda. Para pombe con búsquedas
        SIN mutaciones no es un problema pues estmos hablando de 0.016s

  */

//************************************************************************
//************************************************************************
// ---------------------------- Main code --------------------------------
//************************************************************************
//************************************************************************

var ws=30;           // window size: discrete to real ratio
var processedData;

function main(file, forceReload)
{
    var track="None";
    var nb=5;            // num bins
    var maxSize=400000;  // maximum number of normalized data to store


    destroyAll(false);

    Server.connect(DEBUG_GBV, user, password);

    processedData=uploadFileAndPreprocess(file, forceReload, track, ws, nb, maxSize,2, "False");//stdev and repreprocess the latest 2

    drawing(processedData, ws, maxSize);
}


// ANALYSIS: UPLOAD FILE AND PREPROCESSING
////////////////////////////////////////////
function uploadFileAndPreprocess(file, forceReload, track, ws, nb, maxSize, stdev, repreprocess)
{
    // READ FILE
    //------------
    if(DEBUG_GBV) console.log("\n----- UPLOAD FILE -----");
    if(DEBUG_GBV) console.log("Name of file: "+file.name);

    var startTime=new Date();
    Server.sendFile(file, forceReload);    // passing it to the server side (best solution for >1MB files)
    if (DEBUG_GBV) console.log("Time spent sending: " + (new Date() - startTime) + "ms");


    // PREPROCESSING: NORMALIZE, STATISTICAL DESCRIPTORS, DISCRETIZE
    //----------------------------------------------------------------
    if (DEBUG_GBV) console.log("\n----- PREPROCESSING -----");

    startTime = new Date();
    var processedData = Server.preprocess(file.name, track, ws, nb, maxSize, stdev, repreprocess);
    if (DEBUG_GBV) console.log("Time spent preprocessing: " + (new Date() - startTime) + "ms");

    /* Esto es lo que habría que hacer para invocar al cambio de cromosoma
    startTime = new Date();
    var processedData=Server.getTrack("chromosome2");
    if (DEBUG_GBV) console.log("Time spent getTrack: " + (new Date() - startTime) + "ms");
    */

    return processedData;
}


// DRAWING
////////////////////////////////
function drawing(processedData, ws, maxSize)
{
    if(DEBUG_GBV) console.log("\n----- DRAWING -----");

    var seqServer=processedData.seq;    // just a sampling of about 100K of the original full length sequence
    var fullLength=processedData.fullLength;
    if(DEBUG_GBV) console.log("Length of seqServer:"+seqServer.length+" (full length seq="+fullLength+")");


    // DATALINE 1 (preprocessed data)
    //----------------------------------
    var startTime=new Date();
    dataLine1(seqServer,0,fullLength,fullLength, maxSize, processedData.mean, processedData.stdev, ws, processedData.chromosomes[0]);
    if(DEBUG_GBV) console.log("Time spent dataLine1: "+ (new Date()-startTime)+"ms");

}


// SEARCH POINTS
////////////////////////////////
function searchPoints()
{
    if(globalDL1.drawn)
    {
        var pattern = $('#patternSearch').val();
        var d       = $('#dSearch').val();

        // DATALINE 2: DRAW POINTS
        //----------------------------------
        var numNucleotidesDraw = globalDL1.dim.width;

        //By now working in the first chromosome returned, then we will abilitate a sel box
        console.log("chrom are"+processedData.chromosomes[0]);
        //RODRIGO
        var startTime=new Date();
        var result=Server.search(pattern,d);
        var pp=JSON.parse(result.points[processedData.chromosomes[0]]);
        drawPoints(pp, result.sizePattern, numNucleotidesDraw);
        if(DEBUG_GBV) console.log("Time spent search: "+ (new Date()-startTime)+"ms");

        // GET ENRICHMENT
        var start = new Date().getTime();

        //var annotations=(Server.annotationsGenes("[" + points + "]", "[\"gene\"]", pattern.length * ws, "left", processedData.chromosomes[0]));
        var numMatches=0;

        //This is a much faster code but does not retrieve gene-annotation mappings
        var gis="";
        var annotations=null;
        for(var i=0;i<processedData.chromosomes.length;i++) {
            var points=JSON.parse(result.points[processedData.chromosomes[i]]);
            for( var j=0;j<points.length;j++)
                points[j]*=ws;
            console.log(points.length+" matches in track "+processedData.chromosomes[i]);
            numMatches+=points.length
            gis+=Server.annotationsGenes("[" + points + "]", "[\"gene\"]", result.sizePattern*ws, "left", processedData.chromosomes[i], "True");
            }
         gis=gis.replace(/,$/, "");

/*      //Will solve posterior annotation queries, but it's very time consuming (increases from 1 to 10
        var gis=""//genes of interest
        var annotations={}//for each position, the genes in it
        for(var i=0;i<processedData.chromosomes.length;i++) {
            var points=JSON.parse(result.points[processedData.chromosomes[i]]);
            for( var j=0;j<points.length;j++)
                points[j]*=ws;
            console.log(points.length+" matches in track "+processedData.chromosomes[i]);
            numMatches+=points.length
            //annotations[processedData.chromosomes[i]]=Server.annotationsGenes("[" + points + "]", "[\"gene\"]", result.sizePattern*ws, "left", processedData.chromosomes[i], "False");
            annotations[processedData.chromosomes[i]]=Server.annotationsGenes("[" + points + "]", "[\"any\"]", globalDL2.dim.width, "center", processedData.chromosomes[i], "False");//this here might be too burdening
            gis+=Server.annotationsGenes("[" + points + "]", "[\"gene\"]", result.sizePattern*ws, "left", processedData.chromosomes[i], "True");

        }
        gis=gis.replace(/,$/, "");*/

        console.log("Total number of matches: "+numMatches);
        var end = new Date().getTime();
        setAnnotations(gis,annotations)
        console.log('Getting annotations took: ' + (end - start));


        var start = new Date().getTime();
        var enrichment = Server.enrichment(gis, "fdr", 0.01)
        var end = new Date().getTime();
        console.log('Enrichment analysis took: ' + (end - start));

        drawEnrichment(enrichment, annotations);

        //GET SEQUENCES, MOTIFS, (ALIGNMENT), CONSENSUS
        //NOTE: alignment takes more than 1s if there's >50 sequences! (using the fastest method: kalign)
        var start = new Date().getTime();
        for( var j=0;j<pp.length;j++)
            pp[j]*=ws;
        var response=Server.nucProfile("["+pp+"]",result.sizePattern*ws,globalSeq.track); //TODO: by now only on this chromosome
        var end = new Date().getTime();
        console.log("Sequence analysis took: "+(end-start))
        setSequences(response)
    }
}




// DESTROY ALL
////////////////////////////////
function destroyAll(clear)
{
    // Reset input files
    $("#files").val('');

    // Empty all SVG images
    var images = ["#lineSeq", "#lineSeq2", "#lineSeq3"];

    for(var i=0;i<images.length;i++)
    {
        if($(images[i]).html() != "")
            $(images[i]).empty();
    }

    // Clear console
    if(clear && DEBUG_GBV)
    {
        console.log(new Array(15).join("\n"));
        console.log("Reset all (with File APIs)...");
    }
}


