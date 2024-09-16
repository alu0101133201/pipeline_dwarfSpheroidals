# Small Telescopes Pipeline
##### Sergio Guerra Arencibia
###### Date: 16-09-24

This repository contains the source code of a pipeline implemented for reducing astronomical data from small aperture, large FOV telescopes (amateur telescopes). The purpose of the pipeline is to reduce and produce ultra-deep data from dwarf-spheroidals satellites from the Milky Way.

### Usage

The pipeline expects a configuration file, which contains information of the reduction that is going to be performed (variables, path of the normalisation rings, etc...). 

pipelineName [-h] configurationFile.conf

The pipeline can be used from the pipeline folder itself, although I usually run it from the folder of each of the galaxies. As long as you take care of the paths (path of the normalisation rings, the rootdir of the execution etc...) everything should be fine.

### Folder structure

The structure of the repository is as follows:
.
├── (d) checkScripts 
├── (d) config 
├── (d) pipelineScripts 
├
├── (-) pipeline_LuM_parallel_functions.sh
├── (-) pipeline_LuM_parallel.sh 
├── (-) README.md 
├── (-) template.conf
└── (-) flat_ring_template.txt

* **checkScripts**: Contains python scripts which are not used by the pipeline. They are to be used by the user when checking the data that is going to be reduced (e.g. number of frames, exposure time of frames, airmass-Time plot, etc...)

* **config**: Contains configuration for software used by the pipeline (SExtractor, scamp and swarp mainly). It also contains the astrometry file is created here during the pipeline run. Usually this has not to be touched. 

* **pipelineScripts**: Contains python scripts used by the pipeline. From checking for bad frames to downloading data from decals. This folder is essential.


Then we have the pipeline itself (*pipeline_LuM_parallel.sh*) and the functions which are in another script (*pipeline_LuM_parallel_functions.sh*). Additionally *template.conf* contains the a template with the configuration for the pipeline itself. Here is specified the details of the reduction to be performed (coordinates, type of flat to use, how to model the background, normalisation ring(s), details from the instrument used, etc...). Finally, *flat_ring_template.txt* is a template about how to define the normalisation ring(s).

Take into account that the configuration file (the one corresponding to *template.conf*) is provided to the pipeline as an argument and the normalisation ring(s) are indicated in the configuration file.
