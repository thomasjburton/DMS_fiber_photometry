%% HOUSEKEEPING
%{

This script will iteratively load tdt/synapse data, organise data into
peri-event signals centred on the event, mask out user-defined artefacts,
downsample, fit isosbestic signal to biosensor signal, calculate dFF,
z-score and normalise to a specified baseline window
Parts of this script are adapted from: 
    https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/licking-bout-epoc-filtering/
    https://www.tdt.com/support/matlab-sdk/offline-analysis-examples/fiber-photometry-epoch-averaging-example/
    https://github.com/tjd2002/tjd-shared-code/blob/master/matlab/photometry/FP_normalize.m
Requires MATLAB SDK :
    https://www.tdt.com/docs/sdk/offline-data-analysis/offline-data-matlab/getting-started/


BEFORE RUNNING SCRIPT:
   -create a csv file listing all tdt file/folder names
   -create a csv file listing all epoc labels exactly as they appear in synapse/tdt
   -create csv files specifying the the different stream labels exactly as
    they appear in synapse/tdt, separated by stream type e.g. biosensor
    ("x465A", "x465B") and isosbestic ("x405A", "x405B")
   -create separate csv files detailing the upper and lower thresholds for 
    artefact rejection as well as the temporal cuttoff. This is all done 
    manually at this stage. See README for details


Thomas Burton (t.burton@unsw.edu.au) Jan 2023

%}

%% PREPARING THE WORKSPACE

%specify full path to tdt/synapse files
DATAPATH= fullfile('N:\');

AnimalDataPaths=string(table2cell(readtable('example_filenames.csv', 'ReadVariableNames', false)));   %specify csv file listing all tdt file/folder names
Stream1=string(table2cell(readtable('biosensor_streams.csv', 'ReadVariableNames', false)));         %specify csv file listing each individual stream for the first stream type
Stream2=string(table2cell(readtable('isos_streams.csv', 'ReadVariableNames', false)));              %specify csv file listing each individual stream for the second stream type
Ref_Epocs=string(table2cell(readtable('epoc_labels.csv', 'ReadVariableNames', false)));             %specify csv file listing each individual epoc label

%Parameters for organising, baselining and downsampling data4
TRANGE=[-5, 10];        %extract peri-event data. First value specifies start of window relative to event onset. second value specifies duration of total window
BASELINE = 2.5;         % specify baseline period as number of seconds from the start of epoc 
N = 64;                 % multiplicative for downsampling


%load up details for artefact rejection (user defined, see README)
upper_artefacts=table2array(readtable('upper_artefacts.csv', 'ReadVariableNames', false));
lower_artefacts=table2array(readtable('lower_artefacts.csv', 'ReadVariableNames', false));
cutoff_ends=table2array(readtable('cutoff_ends.csv', 'ReadVariableNames', false));
t_start=8;  %specify number of seconds to discard after the start of a recording. Set to zero if not needed.


%% Pre-processing data to produce peri-event, z-scored and baselined dFF signals 


    
for id=1:length(AnimalDataPaths);       %looping through each data folder in the directory

   
    CurrentID=char(AnimalDataPaths(id));         %variable that can cycle through each data folder in the directory
    BLOCKPATH= fullfile(DATAPATH,CurrentID);   %specify the data directory to be analysed
    
    
    %//////////////////TDT data scripts\\\\\\\\\\\\\\\\\\\\\\\\\
    data = [];    %clear data structure for current iteration
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'epocs', 'scalars', 'streams'}); %run the main function to extract raw data

    for x=1:length(Ref_Epocs)       %looping through each epoc label defined in the tdt file 
        CurrentEpoc=Ref_Epocs(x);   %specify current epoc
        
    %clear working martices
        data_filt = [];
        epocs_working=[];      
        F405cell=[];
        F465cell=[];
        F405cell_filt=[];
        F465cell_filt=[];
        Epoc_filt=[];   
        F405cell_filt_nz=[];
        F465cell_filt_nz=[];
        F405mat = [];

        if sum(strcmp(fieldnames(data.epocs),CurrentEpoc)) > 0      %determine whether current file contains current epoc label
            
        data_filt = TDTfilter(data, CurrentEpoc, 'TIME', TRANGE);   %I have removed the epoc filter so that all events are included
        
        epocs_working(:,1)=data_filt.epocs.(CurrentEpoc).data;      %retrieve raw event data from the data structure and place into a working matrix
        epocs_working(:,2)=data_filt.epocs.(CurrentEpoc).onset;     %as above for event onset
        epocs_working(:,3)=data_filt.epocs.(CurrentEpoc).offset;    %as above for event offset
  

    % store epoc/event data in an array organised by recording (rows) and
    % stream/label within the recording (columns). NOTE: since a pre-event window is specified
    % by TRANGE, any events that occur with a latency less than this window
    % from the start of the recording will not be associated with a
    % photometry stream. To maintain positional equivalency between each
    % event and its corresponding photometry data, the below line removes
    % epocs that will have no peri-event photometry signal calculated due
    % to this issue

        EpocCat_allanimals{id,x}=epocs_working(find(epocs_working(:,3)>=-TRANGE(1,1)),:);

        %specify current stream identity
        dLight = char(Stream1(x));
        Isos = char(Stream2(x));  
 
   
        %extracting peri-event photometry streams 
        F405cell=data_filt.streams.(Isos).filtered;      %extract peri-event data for 405
        F465cell=data_filt.streams.(dLight).filtered;    %as above for 465

        F405min=min(cellfun('length',F405cell));    %calculate minimum raw stream size
        F465min=min(cellfun('length',F465cell));
        F405mat=zeros(length(F405cell),F405min);  %creating matrix of zeros to deal with different sized array elements
        F465mat=zeros(length(F465cell),F465min);


% ////////////////////////// ARTEFACT REMOVAL \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

for i=1:length(F405cell);       %initiate loop that will run through each peri-event cell in array
    
   %clear working martices
    F405mat_filt=[];
    F405mat_hold=[];
    F465mat_hold=[];

    F405mat_hold=cell2mat(F405cell(1,i));               %extract and hold each peri-event 405 signal
    F405mat(i,1:F405min)=F405mat_hold(:,1:F405min);     %store each peri-event 405 signal in a single matrix (only up to min raw stream size)
        
    F465mat_hold=cell2mat(F465cell(1,i));               %as above for 465
    F465mat(i,1:F465min)=F465mat_hold(:,1:F465min);
    
    % identify peri-event windows of the 405 stream that contain artefacts
    % identify any timepoints where the measured signal exceeds or falls short of the manually determined threshold values
    % NOTE: If a given timepoint falls between these values (accepted signal) a 0 will be
    % placed in the corresponding cell. A 1 is inserted if the signal falls outside the bounds 
    % (i.e.artefacts =1). The subsequent conditional statements create a
    % masking matrix (Epoc_filt) that specifies whether an artefact was detected or not
    % for each peri-event 405 stream (col1) or whether the onset of an epoc
    % lies outside the user-defined session bounds (col2)

    F405mat_filt=[];
    F405mat_filt=F405mat>upper_artefacts(id,x) | F405mat<lower_artefacts(id,x); 
    
    
    % Determine whether the current peri-event block for the 405 stream contains a 1 at any point. If so, then mark this peri-event with a 1 (= artefact), 
    % otherwise mark the block with a 0 (= viable block)

   if sum(F405mat_filt(i,:))>0      
       Epoc_filt(i,1)=1;            
   else
       Epoc_filt(i,1)=0;            
   end
   

   %Specifying epocs with onsets that lie outside the user-defined start and end of the
   %recordings here (viable =0; excluded =1). 
   if epocs_working(i,2)+TRANGE(1,1)>t_start & epocs_working(i,2)-TRANGE(1,1)<cutoff_ends(id,x); 
       Epoc_filt(i,2)=0;    
   else
       Epoc_filt(i,2)=1;    
   end
       
   
   if Epoc_filt(i,1)==0 & Epoc_filt(i,2)==0    %determine whether the peri-event block is artefact or viable (in terms of amplitude and time)
       F405cell_filt(i,:)=F405mat(i,1:F405min);   %matrix for viable 405 signals
       F465cell_filt(i,:)=F465mat(i,1:F465min);   %matrix for viable 465 signals
   else
       F405cell_filt(i,1:F405min)=0;   %zeros in matrix for non viable 405 signals
       F465cell_filt(i,1:F465min)=0;   %zeros in matrix for non viable 465 signals
   end
end

    %storing epoc filtering mask and the filtered streams in cell arrays
    Epoc_filt_all{id,x}=Epoc_filt;
    F405cell_filt_all{id,x}=F405cell_filt;          
    F465cell_filt_all{id,x}=F465cell_filt;


% Removing zeros resulting from excising artefacts for both streams
%NOTE: this will create a positional discrepancy between the event/epoc
%stream and the peri-event signal arrays. Only to be used for signal fitting.
for i=1:length(F405cell_filt);
F405cell_filt_nz(:,i)=nonzeros(F405cell_filt(:,i));     
F465cell_filt_nz(:,i)=nonzeros(F465cell_filt(:,i));    
end



% ////////////////////////// DOWNSAMPLING \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

F405cell_filt_nz_ds = zeros(size(F405cell_filt_nz(:,1:N:end-N+1)));
F465cell_filt_nz_ds = zeros(size(F465cell_filt_nz(:,1:N:end-N+1)));
for k = 1:size(F405cell_filt_nz,1);        
    F405cell_filt_nz_ds(k,:) = arrayfun(@(i) mean(F405cell_filt_nz(k,i:i+N-1)),1:N:length(F405cell_filt_nz)-N+1);
    F465cell_filt_nz_ds(k,:) = arrayfun(@(i) mean(F465cell_filt_nz(k,i:i+N-1)),1:N:length(F465cell_filt_nz)-N+1);
end

F405cell_filt_ds = zeros(size(F405cell_filt(:,1:N:end-N+1)));
F465cell_filt_ds = zeros(size(F465cell_filt(:,1:N:end-N+1)));
for k = 1:size(F405cell_filt,1);        
    F405cell_filt_ds(k,:) = arrayfun(@(i) mean(F405cell_filt(k,i:i+N-1)),1:N:length(F405cell_filt)-N+1);
    F465cell_filt_ds(k,:) = arrayfun(@(i) mean(F465cell_filt(k,i:i+N-1)),1:N:length(F465cell_filt)-N+1);
end


% /////////////////// ISOS FITTING AND dFF CALCULATION \\\\\\\\\\\\\\\\\\\\\\\\

% creating working matrices
Y_fit_all_nz=zeros(size(F405cell_filt_nz_ds));
Y_dF_all_nz=zeros(size(F405cell_filt_nz_ds));
dFF_nz=zeros(size(F405cell_filt_nz_ds));

%Fit the 405 onto the 465 stream (polyfit), scale and subtract data to detrend photobleaching
bls = polyfit(F405cell_filt_nz_ds(1:end), F465cell_filt_nz_ds(1:end), 1);
Y_fit_all_nz = bls(1) .* F405cell_filt_nz_ds + bls(2);
Y_dF_all_nz = (F465cell_filt_nz_ds - Y_fit_all_nz);
dFF_nz = (Y_dF_all_nz)./Y_fit_all_nz;


Y_fit_all=zeros(size(F405cell_filt_ds));
Y_dF_all=zeros(size(F405cell_filt_ds));
dFF=zeros(size(F405cell_filt_ds));


%calculating dFF while maintaining raw epoc position 
for i=1:length(F405cell_filt_ds(:,1));
if sum(F405cell_filt_ds(i,:))>0;
Y_fit_all(i,:) = bls(1) .* F405cell_filt_ds(i,:) + bls(2);
Y_dF_all(i,:) = (F465cell_filt_ds(i,:) - Y_fit_all(i,:));
else
    Y_fit_all(i,:) =0;
    Y_dF_all(i,:) =0;
end
end

dFF = (Y_dF_all)./Y_fit_all;


%storing polyfit results and dFF data in cell arrays
bls_all{id,x}=bls;
Y_fit_all_array{id,x}=Y_fit_all;
Y_dF_all_array{id,x}=Y_dF_all;
Y_fit_all_nz_array{id,x}=Y_fit_all_nz;
Y_dF_all_nz_array{id,x}=Y_dF_all_nz;
dFF_all{id,x}=dFF;
dFF_nz_all{id,x}=dFF_nz;


% /////////////////// z-scoring & baselining \\\\\\\\\\\\\\\\\\\\\\\\


    %"local" refers to calculations made for each peri-event signal (i.e Xs
    %either side of event. 

     %creating working matrices
     local_mean_dFF=zeros(length(dFF),1);
     local_std_dFF=zeros(length(dFF),1);
     local_z_dFF_std=zeros(size(dFF));
     
     Local_2_5s_mean_dFF=zeros(length(dFF),1);
     local_2_5s_std_dFF=zeros(length(dFF),1);
     local_2_5s_z_dFF_std=zeros(size(dFF));

     %calculate the mean and stdev for each peri-event signal, then
     %calculate z-score per event
     for i=1:length(dFF(:,1));
        local_mean_dFF(i,1)=mean(dFF(i,:));
        local_std_dFF(i,1)=std(dFF(i,:));
        local_z_dFF_std(i,:)=(dFF(i,:)-local_mean_dFF(i,1))/local_std_dFF(i,1);              
     end

     %store mean, stdev and z-score results in cell arrays
     local_mean_dFF_all{id,x}=local_mean_dFF;
     local_std_dFF_all{id,x}=local_std_dFF;
     local_z_dFF_std_all{id,x}=local_z_dFF_std;


    %specify baselining variables
    bins=length(dFF(1,:));              %specify number of bins/measurements in downsampled stream
    bps=bins/TRANGE(1,2);               %calculate number of bins/measurements per second
    norm_window=round(BASELINE*bps);    %calculate the window for normalising/baselining


     %calculate the mean and stdev for the baselining window for each peri-event signal, then
     %calculate z-score per event using these values
    for i=1:length(dFF(:,1));
        local_2_5s_mean_dFF(i,1)=mean(dFF(i,1:norm_window));    
        local_2_5s_std_dFF(i,1)=std(dFF(i,1:norm_window));
        local_2_5s_z_dFF_std(i,:)=(dFF(i,:)-local_2_5s_mean_dFF(i,1))/local_2_5s_std_dFF(i,1);

    end
    
     %store baseline mean, stdev and resulting z-score results in cell arrays  
     local_2_5s_mean_dFF_all{id,x}=local_2_5s_mean_dFF;
     local_2_5s_std_dFF_all{id,x}=local_2_5s_std_dFF;
     local_2_5s_z_dFF_std_all{id,x}=local_2_5s_z_dFF_std;


        end
        end
        end
        


