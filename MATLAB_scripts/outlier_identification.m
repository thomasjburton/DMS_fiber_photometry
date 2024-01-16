%% HOUSEKEEPING
%{

This script will iteratively load tdt/synapse data and plot raw demodulated streams
Parts of this script are adapted from: 
https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/licking-bout-epoc-filtering/
requires MATLAB SDK https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/getting-started/

BEFORE RUNNING SCRIPT:
  -create a csv file listing all tdt file/folder names
  -create csv files specifying the the different stream labels exactly as
  they appear in synapse/tdt, separated by stream type e.g. biosensor
  ("x465A", "x465B") and isosbestic ("x405A", "x405B")

Thomas Burton (t.burton@unsw.edu.au) Jan 2023

%}

%% PREPARING THE WORKSPACE

%specify full path to tdt/synapse files
DATAPATH= fullfile('N:\');    

AnimalDataPaths=string(table2cell(readtable('example_filenames.csv', 'ReadVariableNames', false)));   %specify csv file listing all tdt file/folder names
Stream1=string(table2cell(readtable('biosensor_streams.csv', 'ReadVariableNames', false)));         %specify csv file listing each individual stream for the first stream type
Stream2=string(table2cell(readtable('isos_streams.csv', 'ReadVariableNames', false)));              %specify csv file listing each individual stream for the second stream type

%provides different colours for each biosensor stream. Isos always plotted in red
colour_pal = [0 0.4 0.8; 0 0.6 0; 1 0.4 0; 0.7 0 0.7; 0.8 0 0];      %1 = blue; 2 = green; 3 = orange; 4 = purple; 5 = red;

%% IMPORTING THE DATA AND PLOTTING


for id=1:length(AnimalDataPaths);       %looping through each data folder in the directory

    CurrentID=char(AnimalDataPaths(id));         %variable that can cycle through each data folder in the directory
    BLOCKPATH= fullfile(DATAPATH,CurrentID);   %specify the data directory to be analysed
    
    
    %//////////////////TDT data scripts\\\\\\\\\\\\\\\\\\\\\\\\\
    data=[];    %clear data structure for current iteration
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'scalars', 'streams'}); %run the main function to extract data
   

for x=1:length(Stream1)    %loop through each channel in the recording

dLight = char(Stream1(x));
Isos = char(Stream2(x));  
 
 %Make a time array based on number of samples (length data stream) and sample freq (fs) of demodulated stream
 time_rawdemod =[];
 time_rawdemod = (1:length(data.streams.(dLight).data))/data.streams.(dLight).fs; 

%plotting raw demodulated data streams

figure('Position',[100, 100, 800, 400])
hold on;
p1 = plot(time_rawdemod, data.streams.(dLight).data,'Color',colour_pal(x,:),'LineWidth',2);
p2 = plot(time_rawdemod, data.streams.(Isos).data,'Color',colour_pal(5,:),'LineWidth',2);
title(strcat({'Raw Demodulated Streams ('},CurrentID,{', '},dLight,{')'}),'fontsize',16);
ylabel('mV','fontsize',16);
axis tight;
legend([p1 p2], {'Biosensor','Isos'});  
end
end





            




