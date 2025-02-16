# Small Telescopes Pipeline: Small FOV, multi-detector adaption
##### Sergio Guerra Arencibia
###### Date: 16-09-24

This repository contains the source code of a pipeline implemented for reducing astronomical data from small aperture and large FOV telescopes. The purpose of the pipeline is to reduce and produce low-surface brightness friendly data.

### Software requirements

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
├ <br />
├── (-) pipeline_LuM_parallel_functions.sh <br />
├── (-) pipeline_LuM_parallel.sh  <br />
├── (-) README.md  <br />
├── (-) template.conf <br />
└── (-) flat_ring_template.txt <br />

* **checkScripts**: Contains python scripts which are not used by the pipeline. They are to be used by the user when checking the data that is going to be reduced (e.g. number of frames, exposure time of frames, airmass-Time plot, etc...)

* **config**: Template for the configuration folder. Contains configuration for software used by the pipeline (SExtractor, scamp and swarp mainly). It also contains the astrometry file is created here during the pipeline run. Usually this has not to be touched. 

* **pipelineScripts**: Contains python scripts used by the pipeline. From checking for bad frames to downloading data from decals. This folder is essential.


Then we have the pipeline itself (*pipeline_LuM_parallel.sh*) and the functions which are in another script (*pipeline_LuM_parallel_functions.sh*). Additionally *template.conf* contains the a template with the configuration for the pipeline itself. Here is specified the details of the reduction to be performed (coordinates, type of flat to use, how to model the background, normalisation ring(s), details from the instrument used, etc...). Finally, *flat_ring_template.txt* is a template about how to define the normalisation ring(s).


##### Notes

Take into account that the configuration file (the one corresponding to *template.conf*) is provided to the pipeline as an argument and the normalisation ring(s) are indicated in the configuration file.

Also, The common normalisation ring (most of the cases will be centered in the image) has to be provided (mandatory) Because it will be used also for selecting what decals bricks are going to be donwloaded for the photometric calibration The 2 rings needed for normalising with them are only requested if the normalisation is going to be done in that way (non mandatory)

## Notes about multi-detector adaptation

This version of the pipeline is an adaptation of the TST and Small telescopes pipeline for multi-ccd arrays such as the WFC of INT or the cameras of LBT. It has only been tested yet in the WFC, but for LBT should be similar. On multi-ccd arrays like Hipercam or Osiris of GTC must be run with caution. The main changes respect to the monolitic pipeline can be summarized as follows:

* **Normalization ring treatement.** It has only been adapted and tested for single, flat ring. It will be adapted for other ring options in the future. The difference with respect to the monolitic lies in the partition of the normalization ring. For that, we have built an script that, given the center of the ring in the mosaic of CCDs, returns different .txt files with the corresponding centers for each one of the CCDs. For WFC of INT, you need to pass the center in the mosaic in X,Y (see documentation of WFC/INT on their webpage). For LBT, since they give an inital astrometry, you can pass the center of the ring in the **FIRST** image in WCS. Then, after astrometrization, a ring must be done again, since sizes change. But here there is a tricky problem. Now, the matrix of CCDs are aligned with the WCS axis. So, if there is a rotated CCD, the file of `pointings_smallGrid` will change NAXIS1 and NAXIS2. Because of that, for those `computeskyForFrame()` with ring that are after `SWARP`, we first check if CCD is rotated with respect to the original .fits or not, and if so we use the `astro-ima` image (which is not rotated) to pass the ring center to RA,DEC and then pass again the RA,DEC into the new system.

* **Multi-layer treatement.** All the pipeline runs so the .fits files created are stored in multi-layer files, one layer for each CCD. In order to keep the important keywords, `AIRMASS` and `DATE_OBS` or `MJD_OBS` are stored in `--hdu=0`.

* **GAIN.** Now the gain might change from one CCD to another. Because of that, we ask for the corresponding keyword in the header instead of for the value.

* **Overscan.** We also add a variable called `$trimsec_key` wich represents the keyword were the illuminated section of the CCD is stored (i.e., the not-overscan region).

* **ASTROMETRY.** Main changes are on the sex-scamp-swarp step. Since they are able to work with multi-layer fits (and scamp works even better), the output catalogs of sextractor are now multi-layer. Swarp configuration changes in order to use a `LANZCOS3` resampling, and keep the temporary files in order to build the pointings folders.

* **Rejection of files.** The current situation is that no rejection is made. This is due to a necessary check in the rejection parameters.

* **Diagnosis and bad files.** The diagnosis and bad files folder is now divided into each of the CCDs, except for the scamp contrast parameters plot.

As with the single detector, multi-detector data should be passed with the `DATA-or` folder structure, but on a multi-layer way where **HDU=0 IS NOT AN IMAGE**.
