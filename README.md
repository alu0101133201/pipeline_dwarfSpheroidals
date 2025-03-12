# Small Telescopes Pipeline
##### Sergio Guerra Arencibia

This repository contains the source code of a pipeline implemented for reducing astronomical data from small aperture and large FOV telescopes. The purpose of the pipeline is to reduce and produce low-surface brightness friendly data.

### Software requirements (this has to be updated)

* gnuastro (currently I'm using 0.22)
* astrometry (using 0.96)
* Scamp (using 2.10)
* Swarp (using 2.41.4)
* SExtractor (using 2.25) 
* funpack (using 1.7.0)

* Astropy (using 5.3.4)

### What do you need to run the pipeline

* Data (raw frames mainly)
* Darks (the pipeline creates a masterdark to correct the BIAS and dark)
* Config directory for software used by the pipeline (template given)
* Configuration file for the specific reduction to perform (.conf file, template give)
* Files specifing the rings to perform the normalisation (also used in data calibration - template given)

##### How the pipeline expects the data

The pipeline is going to look for the following things:

* configuration file (.conf extension) - **Provided to the pipeline when calling it (see *Usage* below)**
    It contains the parameters that belong to the current execution of the pipeline.
    One of the parameters of the file is the ROOTDIR path.

* DATA-or - **Location: path specified in ROOTDIR variable**
    This folder contains the data. The data has to be organised by nights, each of themcalled nightN (with N the number of the night)
    If for whatever reason it doesn't make sense to you to process the data by nights, just place all the data in a *night1* folder

* dark - **Location: path specified in ROOTDIR variable**
    This folder contains the darks. Same formas as the data (i.e. nightN with N the night number)

* Rings file(s) -  **Location: path specified in conf file**
    Definition of ring(s). At least the common ring is mandatory always.
    They are needed for doing two things in the pipeline. The first one is the normalisation of the frames. We use a ring to follow the illumination pattern and avoid the (frequent) central object. The second one is for the photometric calibration. Due to our large FOV we do not download the full decals mosaic, we select some bricks and download them. This bricks are selected based on four different positions of the rings

* config file (configuration) - **Location: path specified in ROOTDIR variable**
    Must contain:
      · Scamp conf file
      · Swarp conf file
      · Sextractor conf files (.conv, .param and .sex)

* filters
    Must contain the transmittances of the filters needed for the reduction.


### Usage

The pipeline expects a configuration file, which contains information of the reduction that is going to be performed (variables, path of the normalisation rings, etc...). 

pipelineName [-h] configurationFile.conf

The pipeline can be used from the pipeline folder itself, although I usually run it from the folder of each of the galaxies. As long as you take care of the paths (path of the normalisation rings, the rootdir of the execution etc...) everything should be fine.

### Folder structure

The structure of the repository is as follows:
.
├── (d) checkScripts <br />
├── (d) config_template <br />
├── (d) pipelineScripts  <br />
├── (d) getFilterCorrection  <br />
├── (d) filters  <br />
├ <br />
├── (-) pipeline_LuM_parallel_functions.sh <br />
├── (-) pipeline_LuM_parallel.sh  <br />
├── (-) README.md  <br />
├── (-) template.conf <br />
└── (-) flat_ring_template.txt <br />

* **checkScripts**: Contains python scripts which are not used by the pipeline. They are to be used by the user when checking the data that is going to be reduced (e.g. number of frames, exposure time of frames, airmass-Time plot, etc...)

* **filters**: Contains the transmittances of the filters needed for the reduction

* **getFilterCorrection**:  Contains python script to be used directly by the user. They allow you to get the colour correction that you need to provide to the pipeline 
in order to take into account the difference in the filters' shape.

* **config**: Template for the configuration folder. Contains configuration for software used by the pipeline (SExtractor, scamp and swarp mainly). It also contains the astrometry file is created here during the pipeline run. Usually this has not to be touched. 

* **pipelineScripts**: Contains python scripts used by the pipeline. From checking for bad frames to downloading data from decals. This folder is essential.


Then we have the pipeline itself (*pipeline_LuM_parallel.sh*) and the functions which are in another script (*pipeline_LuM_parallel_functions.sh*). Additionally *template.conf* contains the a template with the configuration for the pipeline itself. Here is specified the details of the reduction to be performed (coordinates, type of flat to use, how to model the background, normalisation ring(s), details from the instrument used, etc...). Finally, *flat_ring_template.txt* is a template about how to define the normalisation ring(s).


##### Notes

Take into account that the configuration file (the one corresponding to *template.conf*) is provided to the pipeline as an argument and the normalisation ring(s) are indicated in the configuration file.

Also, The common normalisation ring (most of the cases will be centered in the image) has to be provided (mandatory) Because it will be used also for selecting what decals bricks are going to be donwloaded for the photometric calibration The 2 rings needed for normalising with them are only requested if the normalisation is going to be done in that way (non mandatory)

##### IMPORTANT calibration note

Since the calibration factors obtained with PANSTARRS imaging, GAIA spectra and SDDS spectra do NOT completely agree, we have decided to calibrate to GAIA spectra. Thus, we have estimated the aperture needed in PANSTARRS (XRe) to recover magnitudes obtained with GAIA spectra. When doing the tests for estimated this aperture we find that in certain fields we find and offset. For solving this we compute this offset and correct it in each run of the pipeline (thus PANSTARRS always agreeing with GAIA)\\ 
GAIA has been chosen over SDSS because we have more spectra, the calibration is more stable, and we have it in the southern hemisphere. It is true that GAIA sources are quite bright (for TST is fine but would be problematic for other telescopes) but since we only need to calibrate Halpha (much harder to saturate in that band) from bigger telescopes we expect to be fine.\\