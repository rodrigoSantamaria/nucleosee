# Nucleosee

Nucleosome is web-based nucleosome map visual browser backed up by BWT-based pattern searches. It can be used for other kinds of high-throughput genomic data as well.

## Working example
If you want to give it a quick look there is a server set up for tests at http://cpg3.der.usal.es/nucleosee
Try loading samples 972h and Dhta1 and perform some analyses as explained in this video: http://cpg3.der.usal.es/nucleosee/nucleosee.mp4

## Help 
For further help please check http://cpg3.der.usal.es/nucleosee/help.pdf

## Install
If you want to install your own Nucleosee server, it is available as an Nginx Docker container at efialto/nucleosee in the Docker repository.
You basically need to perform 3 steps:
### 1) Install Docker
Visit https://docs.docker.com/install/ to install Docker on your machine.
### 2) Setup host folders
Download the annotations folder at http://vis.usal.es/rodrigo/nucleosee/annotations.zip
Unzip it at your preferred location (ann_path). You can check its folder structure and add 
your own organism annotations.
Optionally, you can download some preprocessed examples at http://vis.usal.es/rodrigo/nucleosee/genomes.zip
Unzip it at your preferred location (gen_path)
This step is optional, as long as you provide two folders with the proper structure for volumes (-v) in the next step, it's fine. More info about it in the help pdf above (section Setting up a server)
### 3) Run Docker container
docker run -it --rm -p 80:80 -v ann_path:/app/annotations -v gen_path:/app/genomes -e SERVERNAME=hostname efialto/nucleosee


## Contact
If you have any doubts about the code or the install, you can write to rodri AT usal DOT es


