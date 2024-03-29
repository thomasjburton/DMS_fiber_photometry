%% HOUSEKEEPING
%{

This script will organise and handle the regression coefficient data
generated in python via the behapy pipeline. csv files that are generated
by Exp2_Reversal_analysis.py can be read into MATLAB and the relevant data
extracted and housed in matrices. Depending on how the design matrix was
built and how groupings were made, the csv will contain subject averaged
peri-event regression coefficients traces/timecourses for fiber photometry
data. For example, if organised by task there will be a single trace per
subject per task for each event type that represents the average regression
coefficient for each timepoint around the event for all such events over
all sessions for the specified task.


BEFORE RUNNING SCRIPT:
   -ensure that there is a csv file generated by Exp2_Reversal_analysis.py
   
Thomas Burton (t.burton@unsw.edu.au) Jan 2023

%}

%% 
%specify csv file to be read
data = readtable('training_beta');

%example code for extracting ipsi and contra lever press events for the RR10 task
ipsi = table2array(data(strcmp(data.Var3, 'ipsilp')&strcmp(data.Var2, 'RR10'),4:end));
contra = table2array(data(strcmp(data.Var3, 'contralp')&strcmp(data.Var2, 'RR10'),4:end));

