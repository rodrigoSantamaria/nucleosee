/*
 ┌────────────────────────────────────────────────────────────┐
 │ upload-file.js                                             │
 ├────────────────────────────────────────────────────────────┤
 │ Description:                                               │
 └────────────────────────────────────────────────────────────┘
 */



// Check if there's support for file reading
if (window.File && window.FileReader && window.FileList && window.Blob)
{
    if(DEBUG_GBV) console.log("Start all (with File APIs)...");

    this.handleEvent=function(evt)
    {
        switch(event.type)
        {
            case 'change':
                var files = evt.target.files;   // FileList object
                var file=files[0];              // By now, just one file

                //NOTE: MD5 computation is just too slow for a large wig (5s for the whole pombe genome)
                //We should replace it by a 'reload' tickmark and just check if the .pic file has been already created otherwise
                /*if(DEBUG_GBV) console.log("\n----- CHECKING FILE -----");
                var startTime=new Date();
                var dfdMd5File = new $.Deferred();
                calculateMD5(dfdMd5File, file);
                dfdMd5File.done(function(hash)
                {
                    if(DEBUG_GBV) console.log("calculateMD5(): "+hash);
                    if(DEBUG_GBV) console.log("Time spent checking (MD5): "+ (new Date()-startTime)+"ms");
                    main(file, hash);
                });*/

                main(file, false) //false if not forcing reload, true otherwise (TODO)

        }
    }
    // Register the custom listener
    document.getElementById('files').addEventListener('change', this, false);
}
else
{
    alert('The File APIs are not fully supported by your browser.');
}




function calculateMD5(dfdMd5File, file)
{
    var reader = new FileReader();
    reader.onloadend = function() {
        dfdMd5File.resolve(md5(reader.result));
    };
    reader.onerror = function() {
        console.error("Could not read the file (MD5)");
    };
    reader.readAsBinaryString(file);
}


