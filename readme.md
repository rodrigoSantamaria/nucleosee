# Nucleosee

Nucleosome is web-based nucleosome map visual browser backed up by BWT-based pattern searches. It can be used for other kinds of high-throughput genomic data as well.

## Working example
If you want to give it a quick look there is a server set up for tests at http://cpg3.der.usal.es/nucleosee
Try loading samples 972h and Dhta1 and perform some analyses as explained in this video: http://cpg3.der.usal.es/nucleosee/nucleosee.mp4

## Help 
For further help please check http://cpg3.der.usal.es/nucleosee/help.pdf

## Install
If you want to install your own Nucleosee server, it is available as an Nginx Docker container at efialto/nucleosee in the Docker repository. This is a stable version ready to run, which may not correspond to the latest code in the GitHub repository, if you want to build a fully updated container you have the components -Dockerfile, etc.- in the repository.

You basically need to perform 3 steps:
### 1) Install Docker
Visit https://docs.docker.com/install/ to install Docker on your machine, if you haven't it already installed.
### 2) Setup host folders
This step is optional, as long as you provide two folders with the proper structure for volumes (-v) in the next step, it's fine. More info about it in the help pdf above (section Setting up a server).

Anyways, we recommend to do it at least the first time you're setting up Nucleosee.

Download the annotations folder at http://vis.usal.es/rodrigo/nucleosee/annotations.zip

Unzip it at your preferred location (`ann_path`). You can check its folder structure and add 
your own organism annotations.

Download some preprocessed examples at http://vis.usal.es/rodrigo/nucleosee/genomes.zip

Unzip it at your preferred location (`gen_path`)


### 3) Run Docker container
```
docker run -it --rm -p 80:80 -v ann_path:/app/annotations -v gen_path:/app/genomes -e SERVERNAME=server_name efialto/nucleosee
```
`gen_path` and `ann_path` relate to the folders created in step 2. `server_name` must be the name of the machine were you're setting up the server (e.g. `signus.opus.uk`)

## Contact
If you have any doubts about the code or the install, you can write to rodri AT usal DOT es

## Citation
Please cite Nucleosee if you find it useful in your research!
  Santamaría, R., Therón, R., Durán, L., García, A., González, S., Sánchez, M., & Antequera, F. (**2019**). Genome-wide search of nucleosome patterns using visual analytics. *Bioinformatics*, 35(13), 2185-2192

