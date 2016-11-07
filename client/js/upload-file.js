/*
 ┌────────────────────────────────────────────────────────────┐
 │ upload-file.js                                             │
 ├────────────────────────────────────────────────────────────┤
 │ Description:                                               │
 └────────────────────────────────────────────────────────────┘
 */



// Check if there's support for file reading...
if (window.File && window.FileReader && window.FileList && window.Blob)
{
    if(DEBUG_GBV) console.log("Start all (with File APIs)...");


    var listenerUploadFile = function(event)
    {
        switch(event.type)
        {
            case 'change':
                var files = event.target.files;   // FileList object
                //var file=files[0];              // By now, just one file
                //checkFile(file);
                checkFile(files);
        }
    };

    // Register the custom listener
    document.getElementById('files').addEventListener('change', listenerUploadFile, false);
}
else
{
    alert('The File APIs are not fully supported by your browser.');
}



/*
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
*/