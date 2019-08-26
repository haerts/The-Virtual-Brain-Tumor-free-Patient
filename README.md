# The-Virtual-Brain-Tumor-free-Patient
Code used for postprocessing analyses in manuscript "Modeling brain dynamics after tumor resection using The Virtual Brain" (Aerts et al. 2019; bioRxiv).

## Workflow

## Paper part 1: Optimizing model parameters after surgery and evaluating stability from pre- to post-operative assessment

(1) Adjust your SC matrices to run with TVBii code (for TVBii code see Schirner et al. 2018 eLife, https://elifesciences.org/articles/28927, https://github.com/BrainModes/The-Hybrid-Virtual-Brain)
> Generate_TVBii_Input_v4_post
  
(2) Run TVBii code
> TVBii_G-Ji_prepnrun_t2.sh

(3) Perform parameter space exploration to identify optimal model parameters
> PSE.m

(4) Postprocessing:
> GTA_SC.m: compute structural network topology measures
> postpro_prep_20190821.R: descriptive analyses
> postpro_regression_20190821.R: assess associations between fitted model parameters, structural network topology measures and cognitive performance measures


## Paper part 2: Virtual neurosurgery proof of concept analyses

(1) Adjust your SC matrices to run with TVBii code
> Generate_TVBii_Input_v4_SS3Tpre: to be used on SC matrix reconstruced using SS3T
> Generate_TVBii_Input_v4_SS3TpreVS: to be used on SC matrix reconstruced using SS3T, after virtual neurosurgery (ie removal of white matter tracts intersecting resection mask)

(2) Run TVBii code: similar to code used for part 1

(3) Perform parameter space exploration to identify optimal model parameters
> PSE_SS3Tpre.m

(4) Evaluation of computational models for prediction of post-surgical brain dynamics
> VSevaluation.m: evaluation
> ModelFit_20190823.R: visualization


