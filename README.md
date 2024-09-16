# Small Telescopes Pipeline
##### Sergio Guerra Arencibia
###### Date: 16-09-24

This repository contains the source code of a pipeline implemented for reducing astronomical data from small aperture, large FOV telescopes (amateur telescopes). The purpose of the pipeline is to reduce and produce ultra-deep data from dwarf-spheroidals satellites from the Milky Way.

### Usage

pipelineName [-h] configurationFile.conf

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
└── (-) template.conf

* **checkScripts**: Contains python scripts which are not used by the pipeline. They are to be used by the user when checking the data that is going to be reduced (e.g. number of frames, exposure time of frames, airmass-Time plot, etc...)

* **config**: Contains configuration for software used by the pipeline (SExtractor, scamp and swarp mainly). It also contains the definition of the normalisation ring and an astrometry file is created here during the pipeline run. Usually this has not to be touched. 

* **pipelineScripts**: Contains python scripts used by the pipeline. From checking for bad frames to downloading data from decals. This folder is essential.


Then we have the pipeline itself (*pipeline_LuM_parallel.sh*) and the functions which are in another script (*pipeline_LuM_parallel_functions.sh*).

Finally *template.conf* contains the a template with the configuration for the pipeline itself. Here is specified the details of the reduction to be performed (coordinates, type of flat to use, how to model the background, details from the instrument used, etc...). 
