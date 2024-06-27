# DMS_fiber_photometry

This folder contains files that were used to process and analyse fiber photometry data detailed in [this paper](https://doi.org/10.1016/j.celrep.2024.113828).

Data accessible [here](https://doi.org/10.6084/m9.figshare.19083647.v2).

This is a transitory approach that allowed us to wrangle, process and analyse photometry data but our team has since developed another framework for handling such data. These MATLAB scripts will therefore not be developed or improved. Going forward, we will be using the [behapy](https://github.com/crnolan/behapy) toolkit. 

The data was generated using a TDT system (RZ5P Processor, Synapse, Tucker-Davis Technologies, TDT) and these scripts use functions provided by MATLAB SDK. We simultaneously recorded two streams of data per recording site: fluorescence emissions from 465 nm (biosensor) and 405 nm (isosbestic) light stimulation. Anywhere from 1-4 subjects were recorded simultaneously (in operant boxes A-D), with data from all subjects per recording packaged into the same raw files. The MATLAB workspace for handling the raw data is organised such that each row is a recording session and each column is a subject within that recording session (also referred to subsequently as a single recording). Subjects were run in operant boxes A-D and the data generated were organised in the MATLAB workspace as columns 1-4. Events (timestamped via a TTL pulse from med-pc) and the two photometry streams for each individual recording site were lablled thus:

Operant Box A (col 1): PrtA ("epocs" or events), x405A (isosbestic stream), x465A (biosensor stream)

Operant Box B (col 2): PrtB ("epocs" or events), x405B (isosbestic stream), x465B (biosensor stream)

etc...

## Files

*outlier_identification.m*
<br />
This script can be run to batch extract photometry data from tdt files and plot raw signals (405 and 465 pairs). The user must specify the filenames to be handled (see example file *example_filenames.csv*) and the order these are specified will determine the order in which the data is positioned in rows. The user can then manually determine the upper and lower thresholds for artefact removal for each recording (we did this for the 405 stream alone). A cutoff for the end of the session can also be determined manually. *outlier_detection_example.jpg* illustrates manually defined upper (8) and lower (5) thresholds and a cutoff for session end (1938.86s). These are entered into separate csv files that reflect the structure of MATLAB workspace (i.e. recordings from boxes A-D represented in cols 1-4). Example files are provided (*cutoff_ends.csv*, *lower_artefacts.csv*, *upper_srtefacts.csv*). The data in (2,3) corresponds to Operant Box C (col3) in the second recording (row2; specified by the order of filenames). 

*preprocessing.m*
<br />
Once the artefact thresholds have been manually determined (and csv files created) for each recording, the data can be batch pre-processed. This script segments the data into user-specified windows centred around events/epocs (e.g. 5s before and 5s after a lever press). Each timepoint in each peri-event trace is assessed for artefact violations. If the value of the 405 at any point in the trace lies outside the user-specified thresholds, the entire trace is discarded. Data is downsampled (with a user specified multiplicative) and the isosbestic 405 is regressed onto the 465 stream using a least-squares linear fit model (*polyfit*) for each recording session. Data is z-scored for each peri-event trace using the mean and stdev of that particular trace within a user-specified baseline window (for each timepoint in the trace: [dFF-mean_dFF]/std_dFF). In this case, it was from -5:-2.5s before each event. Variables and data of interest are stored in cell arrays with their position reflecting recording order and Operant box.

*event_handling.m*
<br />
This script requires a specification of the identities of the coded events. It is specific to the experimental situation characterising this study. The data provided on FigShare is packaged into a single MATLAB structure for each experiment, with everything now ordered by subject and by session. This script provides code for extracting this data. This script will iteratively cycle through each recording and organise the photometry data (z-scored and baselined dFF) into different basic event types. 

*event_sequences.m*
<br />
This script requires a specification of the identities of the coded events. It is specific to the experimental situation characterising this study. The data provided on FigShare is packaged into a single MATLAB structure for each experiment, with everything now ordered by subject and by session. This script is designed to follow the initial event_handling, and provides code for extracting events categorized accoring to the events that preceded them (being either another press, a magazine entry or neither of these events for a pre-determined period of time). This script will iteratively cycle through each recording and organise the photometry data (z scored and baselined dFF) into different event sequence categories. At the end of this script is a small example section showing how the area under the curve (AUC) was calculated.

*Learning_Simulations.m*
<br />
This script is designed to follow the initial *event_handling.m*, and uses the pre-extracted events to create and simulate a learning rule that updates across and within sessions as animals make lever presses and magazine entries. This script will iteratively cycle through each recording session for each animal in order (from first training session to last).

*REPOSITORY_CDeg_code.m*
<br />
This script is a modification of the *event_handling.m* script, which will will distinguish between free and earned outcome deliveries for contingency degradation. This script will iteratively cycle through each recording and organise the photometry data (z-scored and baselined dFF) into different event types. 

*regcoeff_readcsv.m*
<br />
This script will read csv files generated by ```regression_preprocessing/Exp2_Reversal_analysis.py```. Subject-based regression coefficient data can be extracted (example provided) for further handling in MATLAB.

<br />
<br />
Parts of these scripts are adapted from: 

    https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/licking-bout-epoc-filtering/
    
    https://www.tdt.com/support/matlab-sdk/offline-analysis-examples/fiber-photometry-epoch-averaging-example/
    
    https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m
    
Requires MATLAB SDK :

    https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/getting-started/
<br />


Waveform analysis was achieved using philjrdb's [ERTsimulation](https://github.com/philjrdb/ERTsimulation).
<br />
<br />
Thomas Burton, Jan 2023
