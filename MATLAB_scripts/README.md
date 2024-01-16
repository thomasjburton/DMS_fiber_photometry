# DMS_fiber_photometry

This folder contains files that were used to process and analyse fiber photometry data detailed in this paper: https://doi.org/10.1101/2022.01.31.478585

This is a transitory approach that allowed us to wrangle, process and analyse photometry data but our team has since developed another framework for handling such data. These MATLAB scripts will therefore not be developed or improved. Going forward, we will be using the following toolkit: https://github.com/crnolan/behapy

The data was generated using a TDT system (RZ5P Processor, Synapse, Tucker-Davis Technologies, TDT) and these scripts use functions provided by MATLAB SDK. We simultaneously recorded two streams of data per recording site: fluorescence emissions from 465 nm (biosensor) and 405 nm (isosbestic) light stimulation. Anywahere from 1-4 subjects were recorded at a given time (operant boxes A-D), with data from all subjects per recording packaged in a single raw file. The MATLAB workspace is organised such that each row is a recording session and each column is a subject within that recording session (also referred to as a single recording). Subjects were run in operant boxes A-D, organised in the MATLAB workspace as columns 1-4. Events (timestamped via a TTL pulse from med-pc) and the two photometry streams for each individual recording site were lablled according to the following:

Operant Box A (col 1): PrtA ("epocs" or events), x405A (isosbestic stream), x465A (biosensor stream)
Operant Box B (col 2): PrtB ("epocs" or events), x405B (isosbestic stream), x465B (biosensor stream)
etc...


*outlier_identification.m*
<br />
This script can be run to batch extract photometry data from tdt files and plot raw signals. The user can then manually determine the upper and lower thresholds for artefact removal for each recording. A cutoff for the end of the session can also be determined manually. *outlier_detection_example.jpg* illustrates manually defined upper (8) and lower (5) thresholds and a cutoff for session end (1938.86s). These are entered into separate csv files that reflect the structure of MATLAB workspace (i.e. recordings from boxes A-D represented in cols 1-4). Example files are provided.





Parts of this script are adapted from: 
    https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/licking-bout-epoc-filtering/
    https://www.tdt.com/support/matlab-sdk/offline-analysis-examples/fiber-photometry-epoch-averaging-example/
    https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m
Requires MATLAB SDK :
    https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/getting-started/


Thomas Burton (t.burton@unsw.edu.au) Jan 2023
