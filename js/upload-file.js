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

                if(DEBUG_GBV) console.log("\n----- CHECKING FILE -----");
                var startTime=new Date();
                if(false)
                {
                    var dfdMd5File = new $.Deferred();
                    calculateMD5(dfdMd5File, file);
                    dfdMd5File.done(function (hash) {
                        if (DEBUG_GBV) console.log("calculateMD5(): " + hash);
                        if (DEBUG_GBV) console.log("Time spent checking (MD5): " + (new Date() - startTime) + "ms");
                        main(file, hash);
                    });
                }
                else
                {
                    var hash="WITHOUT_MD5";
                    if (DEBUG_GBV) console.log("calculateMD5(): " + hash);
                    if (DEBUG_GBV) console.log("Time spent checking (MD5): " + (new Date() - startTime) + "ms");
                    main(file, hash);
                }

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


