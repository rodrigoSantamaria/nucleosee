

// Check if there's support for file reading
if (window.File && window.FileReader && window.FileList && window.Blob)
    {
    console.log("Start all (with File APIs)...");
    this.handleEvent=function(evt)
    {
        switch(event.type)
            {
            case 'change':
                var files = evt.target.files; // FileList object
                file=files[0]  //By now, just one file

                main(file);
            }
    }
    // Register the custom listener
    document.getElementById('files').addEventListener('change', this, false);
    }
else
    {
    alert('The File APIs are not fully supported by your browser.');
    }
