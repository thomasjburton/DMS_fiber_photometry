
%% HOUSEKEEPING
%{

This script will organise pre-processed fiber photometry data (organised by subject and session) 
into event types. This script was written for analysing a free-operant lever-pressing task using
med-associates/med-pc gear. Rats were able to make 2 lever press responses
(left and right; LLP and RLP) that were situated either side of a recessed reward
magazine. Lever-pressing resulted in the delivery of distinct reward
outcomes (grain pellets or sucrose solution; PEL and SUC). Onsets and offsets for each
lever press response, magazine entry (via IR beam break) and outcome delivery were recorded and each event
occurrence has a corresponding stream of data centred on the event.

BEFORE RUNNING SCRIPT:
   -create a data structure containing positionally mapped event and pre-processed
    photometry data.
   


Thomas Burton (t.burton@unsw.edu.au) Jan 2023

%}


%% Provide the epoc-event identity mapping

LLP=1;          %left lever press
RLP=4;          %right lever press
PEL=16;         %pellet delivery
SUC=8;          %sucrose delivery
ENTRY_CODE=2;   %magazine entry

%% Extracting data from archiving structure

%using example data structure "Exp2_Reversal"
load('Exp2_Reversal.mat')

for i=1:length(Exp2_Reversal)
events(i,:)=Exp2_Reversal(i).events;
events_mask(i,:)=Exp2_Reversal(i).events_mask;
signals(i,:)=Exp2_Reversal(i).z_base_dFF;
end


%% 


for h=1:length(events(:,1))
    
    for x=1:length(events(1,:))

    %clear working matrices
    events_working=[];
    events_filt_working=[];
    Signals_working_dFF=[];

   
    %retrieve current event and signal data from cell array
    events_working=events{h,x};  
    events_filt_working=events_mask{h,x};
    Signals_working_dFF=signals{h,x};    

    if length(events_working) >0
    %adjusting the med-pc event codes output to reflect event onset and
    %offset (col5). Also including event number (col2)
        events_working(1,5)=events_working(1,1);    %manually inserting for first cell only
        events_working(1,2)=1;                      %manually inserting for first cell only
        for i=2:length(events_working);
        events_working(i,5)=events_working(i,1)-events_working(i-1,1);
        events_working(i,2)=i;
        end
    

%////////////////////////////////////////////////////////////////////    
%//////////////Broad Behavioural Event Categories\\\\\\\\\\\\\\\\\\\\
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

   %clear working matrices    
    Locate_LLP=[];
    Locate_RLP=[];
    Locate_Mag=[];
    Locate_Pel=[];
    Locate_Suc=[];
    Locate_Mag=[];
    Locate_LLP_off=[];
    Locate_RLP_off=[];
    Locate_Mag_off=[];
    LLP_Signal_dFF=[];
    RLP_Signal_dFF=[];
    Mag_Signal_dFF=[];
    First_Mag_Pel=[];
    First_Mag_Suc=[];
    Mag_working=[];
    First_Mag_Pel_Signal_dFF=[];
    First_Mag_Suc_Signal_dFF=[];
    Locate_Mag_Empty=[];

    %Locate event numbers for each event type
    Locate_LLP=find(events_working(:,5)==LLP) ;  
    Locate_RLP=find(events_working(:,5)==RLP) ;
    Locate_Mag=find(events_working(:,5)==ENTRY_CODE) ;
    Locate_Pel=find(events_working(:,5)==PEL) ;
    Locate_Suc=find(events_working(:,5)==SUC) ;
    Locate_LLP_off=find(events_working(:,5)==-LLP) ;
    Locate_RLP_off=find(events_working(:,5)==-RLP) ;
    Locate_Mag_off=find(events_working(:,5)==-ENTRY_CODE) ;
   
    %store event numbers for each event type in cell arrays
    Locate_LLP_all{h,x}=Locate_LLP;
    Locate_RLP_all{h,x}=Locate_RLP;
    Locate_Mag_all{h,x}=Locate_Mag;
    Locate_Pel_all{h,x}=Locate_Pel;
    Locate_Suc_all{h,x}=Locate_Suc;
    Locate_LLP_off_all{h,x}=Locate_LLP_off;
    Locate_RLP_off_all{h,x}=Locate_RLP_off;
    Locate_Mag_off_all{h,x}=Locate_Mag_off;
   


%extraction of signal data under broad categories of behavioural data
for k=1:length(Locate_LLP);
    if Locate_LLP(k,1)<= length(Signals_working_dFF(:,1)); 
   LLP_Signal_dFF(k,:)=Signals_working_dFF(Locate_LLP(k,1),:);    
    end
end

for k=1:length(Locate_RLP);
   if Locate_RLP(k,1)<= length(Signals_working_dFF(:,1)); 
     RLP_Signal_dFF(k,:)=Signals_working_dFF(Locate_RLP(k,1),:);    
   end
end
    
for k=1:length(Locate_Mag);
   if Locate_Mag(k,1)<= length(Signals_working_dFF(:,1)); 
     Mag_Signal_dFF(k,:)=Signals_working_dFF(Locate_Mag(k,1),:);    
   end
end

for k=1:length(Locate_Pel);
   if Locate_Pel(k,1)<= length(Signals_working_dFF(:,1)); 
     Pel_Signal_dFF(k,:)=Signals_working_dFF(Locate_Pel(k,1),:);    
   end
end

for k=1:length(Locate_Suc);
   if Locate_Suc(k,1)<= length(Signals_working_dFF(:,1)); 
     Suc_Signal_dFF(k,:)=Signals_working_dFF(Locate_Suc(k,1),:);    
   end
end

%Store broad event signals in cell arrays
LLP_Signal_all{h,x}=LLP_Signal_dFF;
RLP_Signal_all{h,x}=RLP_Signal_dFF;
Mag_Signal_all{h,x}=Mag_Signal_dFF;
Pel_Signal_all{h,x}=Pel_Signal_dFF;
Suc_Signal_all{h,x}=Suc_Signal_dFF;




%create a working variable that removes all event data except that pertaining 
%to mag entry onset.
    
    for i=1:length(events_working);
    if events_working(i,5) == 2;
        Mag_working(i,1:5)=events_working(i,1:5);
    else Mag_working(i,1:5)=0;
    end 
    end


%Locating Mag entries immediately following outcome delivery (i.e. contact with reward)
   for j=1:size(Locate_Pel)-1;
       if Locate_Pel(j)<Locate_Mag(end) & Locate_Pel(j+1)>find(Mag_working(:,2)>Locate_Pel(j),1);
    First_Mag_Pel(j,1)=find(Mag_working(:,2)>Locate_Pel(j),1);
       end 
   end

   if length(Locate_Pel)>0;
   j=length(Locate_Pel);
   if Locate_Pel(j)<Locate_Mag(end)
    First_Mag_Pel(j,1)=find(Mag_working(:,2)>Locate_Pel(j),1);
   end
  First_Mag_Pel=nonzeros(First_Mag_Pel);
   end
   
   for j=1:size(Locate_Suc)-1;
       if Locate_Suc(j)<Locate_Mag(end) & Locate_Suc(j+1)>find(Mag_working(:,2)>Locate_Suc(j),1);
    First_Mag_Suc(j,1)=find(Mag_working(:,2)>Locate_Suc(j),1);
       end 
   end

   if length(Locate_Suc)>0;
   j=length(Locate_Suc);
   if Locate_Suc(j)<Locate_Mag(end)
    First_Mag_Suc(j,1)=find(Mag_working(:,2)>Locate_Suc(j),1);
   end
   First_Mag_Suc=nonzeros(First_Mag_Suc);
   end
  
 %store event numbers for reward contact
 First_Mag_Pel_all{h,x}=First_Mag_Pel;
 First_Mag_Suc_all{h,x}=First_Mag_Suc;



%extract signals around reward contact
for k=1:length(First_Mag_Pel);
   if First_Mag_Pel(k,1)<= length(Signals_working_dFF(:,1)); 
     First_Mag_Pel_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Pel(k,1),:);    
   end
end

for k=1:length(First_Mag_Suc);
   if First_Mag_Suc(k,1)<= length(Signals_working_dFF(:,1)); 
     First_Mag_Suc_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Suc(k,1),:);    
   end
end

%store signal data for reward contact
First_Mag_Pel_Signal_all{h,x}=First_Mag_Pel_Signal_dFF;
First_Mag_Suc_Signal_all{h,x}=First_Mag_Suc_Signal_dFF;



%locating empty mag entries (i.e. unrewarded)
Locate_Mag_Empty=setdiff(Locate_Mag,First_Mag_Pel);
Locate_Mag_Empty=setdiff(Locate_Mag_Empty,First_Mag_Suc);

%store empty mag events in cell array
Locate_Mag_Empty_all{h,x}=Locate_Mag_Empty;


%/////////////////////////////////////////////////////////////////////    
%///////////////////////Sequences of Behaviour\\\\\\\\\\\\\\\\\\\\\\\\
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


%clear working matrices
LLP_LLP=[];         %two consecutive left lever presses
LLP_Mag=[];         %left lever press followed by all mag entries
LLP_Mag_Empty=[];   %left lever press followed by an unrewarded mag entry
LLP_Mag_Pel=[];     %left lever press followed by a mag entry with a pel reward
LLP_Mag_Suc=[];     %etc, etc, etc.
RLP_RLP=[];
RLP_Mag=[];
RLP_Mag_Empty=[];
RLP_Mag_Pel=[];
RLP_Mag_Suc=[];
LLP_RLP=[];
RLP_LLP=[];

%//////////////////////////////LLPs////////////////////////////////////////

for j=1:size(Locate_LLP)-1; %not looking at final LP of session
    %find LLPs that occur before the final mag entry of the session and where 
    %the subsequent LLP occurs prior to a mag entry
       if Locate_LLP(j)<Locate_Mag(end) & Locate_LLP(j+1)<find(Mag_working(:,2)>Locate_LLP(j),1); 
    LLP_LLP(j,1)=Locate_LLP(j);
    %added this code in to ensure that instances of LL-RL_Mag (as might occur at the end of a training block) 
       elseif Locate_LLP(j)<Locate_Mag(end) & Locate_RLP(find(Locate_RLP(:,1)>Locate_LLP(j),1))<find(Mag_working(:,2)>Locate_LLP(j),1) & Locate_LLP(j+1)>Locate_RLP(find(Locate_RLP(:,1)>Locate_LLP(j),1))
       LLP_RLP(j,1)=Locate_LLP(j);
       end 
end

LLP_LLP=nonzeros(LLP_LLP);
LLP_RLP=nonzeros(LLP_RLP);
LLP_Mag=setdiff(Locate_LLP,LLP_LLP);
LLP_Mag=setdiff(LLP_Mag,LLP_RLP);

LLP_LLP_all{h,x}=LLP_LLP;
LLP_Mag_all{h,x}=LLP_Mag;


for j=1:size(LLP_Mag); %not looking at final LP of session
    %create a working variable that will locate the first occurence of each mag entry type for a give LP-Mag combination
    %note that this covers all combinations of LP and outcome so one doesn't have to specify which press earned which outcome. 
    M_O_working=zeros(3,2); 
    
    if length(Locate_Mag_Empty)>0;
    if LLP_Mag(j)<Locate_Mag_Empty(end);
    M_O_working(1,1)=find(Locate_Mag_Empty(:,1)>LLP_Mag(j),1); %locate the first empty mag entry directly after an LLP. Record position in Locate_Mag_Empty vector
    M_O_working(1,2)=Locate_Mag_Empty(M_O_working(1,1),1);      %record the event cell position reference number for this mag entry
    end
    end
    
    %as above for mag entries containing sucrose
    if length(First_Mag_Suc)>0;
    if LLP_Mag(j)<First_Mag_Suc(end);
    M_O_working(2,1)=find(First_Mag_Suc(:,1)>LLP_Mag(j),1);
    M_O_working(2,2)=First_Mag_Suc(M_O_working(2,1),1);
    end
    end
    
    %as above for mag entries containing a pellet
    if length(First_Mag_Pel)>0;
    if LLP_Mag(j)<First_Mag_Pel(end);
    M_O_working(3,1)=find(First_Mag_Pel(:,1)>LLP_Mag(j),1);
    M_O_working(3,2)=First_Mag_Pel(M_O_working(3,1),1);
    end
    end
    
    %determine which mag entry type immediately followed the LLP. The
    %lowest value in the second column will reveal which mag entry type was
    %closest to the LLP followed by mag entry
    if M_O_working(1,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(1,2)>0;
    LLP_Mag_Empty(j,1)=LLP_Mag(j,1);
    elseif M_O_working(2,2)== min(nonzeros(M_O_working(:,2))) & M_O_working(2,2)>0;
    LLP_Mag_Suc(j,1)=LLP_Mag(j,1);
    elseif M_O_working(3,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(3,2)>0;
    LLP_Mag_Pel(j,1)=LLP_Mag(j,1);
    end
end


LLP_Mag_Empty=nonzeros(LLP_Mag_Empty);
LLP_Mag_Pel=nonzeros(LLP_Mag_Pel);
LLP_Mag_Suc=nonzeros(LLP_Mag_Suc);


%finding the empty mag entries following LPs to plot data aligned to this
%point rather than the LP
magaligned_working=[];
magaligned_LLP_Mag_Empty=[];
for j=1:length(LLP_Mag_Empty);
    magaligned_working(j,1)=find(Locate_Mag_Empty(:,1)>LLP_Mag_Empty(j),1);
    magaligned_LLP_Mag_Empty(j,1)=Locate_Mag_Empty(magaligned_working(j,1),1);
end


LLP_Mag_Empty_all{h,x}=LLP_Mag_Empty;
LLP_Mag_Pel_all{h,x}=LLP_Mag_Pel;
LLP_Mag_Suc_all{h,x}=LLP_Mag_Suc;
magaligned_LLP_Mag_Empty_all{h,x}=magaligned_LLP_Mag_Empty;



%//////////////////////////////RLPs////////////////////////////////////////
%as above for RLPs
for j=1:size(Locate_RLP)-1; %not looking at final LP of session
       if Locate_RLP(j)<Locate_Mag(end) & Locate_RLP(j+1)<find(Mag_working(:,2)>Locate_RLP(j),1);
    RLP_RLP(j,1)=Locate_RLP(j);
    elseif Locate_RLP(j)<Locate_Mag(end) & Locate_LLP(find(Locate_LLP(:,1)>Locate_RLP(j),1))<find(Mag_working(:,2)>Locate_RLP(j),1) & Locate_RLP(j+1)>Locate_LLP(find(Locate_LLP(:,1)>Locate_RLP(j),1))
       RLP_LLP(j,1)=Locate_RLP(j);
         end 
end
  
RLP_RLP=nonzeros(RLP_RLP);
RLP_LLP=nonzeros(RLP_LLP);
RLP_Mag=setdiff(Locate_RLP,RLP_RLP);
RLP_Mag=setdiff(RLP_Mag,RLP_LLP);

RLP_RLP=nonzeros(RLP_RLP);
RLP_Mag=setdiff(Locate_RLP,RLP_RLP);

RLP_RLP_all{h,x}=RLP_RLP;
RLP_Mag_all{h,x}=RLP_Mag;

for j=1:size(RLP_Mag); %not looking at final LP of session
    M_O_working=zeros(3,2); %create a working variable that will locate the first occurence of a mag entry type for a give LP-Mag combination
    
    if length(Locate_Mag_Empty)>0;
    if RLP_Mag(j)<Locate_Mag_Empty(end);
    M_O_working(1,1)=find(Locate_Mag_Empty(:,1)>RLP_Mag(j),1);
    M_O_working(1,2)=Locate_Mag_Empty(M_O_working(1,1),1);
    end
     end
     
     if length(First_Mag_Suc)>0;
    if RLP_Mag(j)<First_Mag_Suc(end);
    M_O_working(2,1)=find(First_Mag_Suc(:,1)>RLP_Mag(j),1);
    M_O_working(2,2)=First_Mag_Suc(M_O_working(2,1),1);
    end
     end
   
      if length(First_Mag_Pel)>0;
    if RLP_Mag(j)<First_Mag_Pel(end);
    M_O_working(3,1)=find(First_Mag_Pel(:,1)>RLP_Mag(j),1);
    M_O_working(3,2)=First_Mag_Pel(M_O_working(3,1),1);
    end
      end
    
    if M_O_working(1,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(1,2)>0;
    RLP_Mag_Empty(j,1)=RLP_Mag(j,1);
    elseif M_O_working(2,2)== min(nonzeros(M_O_working(:,2))) & M_O_working(2,2)>0;
    RLP_Mag_Suc(j,1)=RLP_Mag(j,1);
    elseif M_O_working(3,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(3,2)>0;
    RLP_Mag_Pel(j,1)=RLP_Mag(j,1);
    end
end
     
RLP_Mag_Empty=nonzeros(RLP_Mag_Empty);
RLP_Mag_Pel=nonzeros(RLP_Mag_Pel);
RLP_Mag_Suc=nonzeros(RLP_Mag_Suc);


magaligned_working=[];
magaligned_RLP_Mag_Empty=[];
for j=1:length(RLP_Mag_Empty);
    magaligned_working(j,1)=find(Locate_Mag_Empty(:,1)>RLP_Mag_Empty(j),1);
    magaligned_RLP_Mag_Empty(j,1)=Locate_Mag_Empty(magaligned_working(j,1),1);
end

RLP_Mag_Empty_all{h,x}=RLP_Mag_Empty;
RLP_Mag_Pel_all{h,x}=RLP_Mag_Pel;
RLP_Mag_Suc_all{h,x}=RLP_Mag_Suc;
magaligned_RLP_Mag_Empty_all{h,x}=magaligned_RLP_Mag_Empty;



%/////////////////Extracting signals///////////////////////////////////////

%clear working matrices
LLP_LLP_Signal_dFF=[];
LLP_Mag_Empty_Signal_dFF=[];
LLP_Mag_Pel_Signal_dFF=[];
LLP_Mag_Suc_Signal_dFF=[];
magaligned_LLP_Mag_Pel_Signal_dFF=[];
magaligned_LLP_Mag_Suc_Signal_dFF=[];
magaligned_LLP_Mag_Empty_Signal_dFF=[];

RLP_RLP_Signal_dFF=[];
RLP_Mag_Empty_Signal_dFF=[];
RLP_Mag_Pel_Signal_dFF=[];
RLP_Mag_Suc_Signal_dFF=[];
magaligned_RLP_Mag_Pel_Signal_dFF=[];
magaligned_RLP_Mag_Suc_Signal_dFF=[];
magaligned_RLP_Mag_Empty_Signal_dFF=[];


%Aligning press-entry sequence signal data to lever press 
for k=1:length(LLP_LLP);
    if LLP_LLP(k,1)<= length(Signals_working_dFF(:,1));
   LLP_LLP_Signal_dFF(k,:)=Signals_working_dFF(LLP_LLP(k,1),:);    
    end
end

for k=1:length(LLP_Mag_Empty);
    if LLP_Mag_Empty(k,1)<= length(Signals_working_dFF(:,1));
   LLP_Mag_Empty_Signal_dFF(k,:)=Signals_working_dFF(LLP_Mag_Empty(k,1),:);    
    end
end

for k=1:length(LLP_Mag_Pel);
    if LLP_Mag_Pel(k,1)<= length(Signals_working_dFF(:,1));
   LLP_Mag_Pel_Signal_dFF(k,:)=Signals_working_dFF(LLP_Mag_Pel(k,1),:);    
    end
end

for k=1:length(LLP_Mag_Suc);
    if LLP_Mag_Suc(k,1)<= length(Signals_working_dFF(:,1));
   LLP_Mag_Suc_Signal_dFF(k,:)=Signals_working_dFF(LLP_Mag_Suc(k,1),:);    
    end
end


for k=1:length(RLP_RLP);
    if RLP_RLP(k,1)<= length(Signals_working_dFF(:,1)); 
   RLP_RLP_Signal_dFF(k,:)=Signals_working_dFF(RLP_RLP(k,1),:);    
    end
end

for k=1:length(RLP_Mag_Empty);
    if RLP_Mag_Empty(k,1)<= length(Signals_working_dFF(:,1)); 
   RLP_Mag_Empty_Signal_dFF(k,:)=Signals_working_dFF(RLP_Mag_Empty(k,1),:);    
    end
end

for k=1:length(RLP_Mag_Pel);
    if RLP_Mag_Pel(k,1)<= length(Signals_working_dFF(:,1));  
   RLP_Mag_Pel_Signal_dFF(k,:)=Signals_working_dFF(RLP_Mag_Pel(k,1),:);    
    end
end

for k=1:length(RLP_Mag_Suc);
    if RLP_Mag_Suc(k,1)<= length(Signals_working_dFF(:,1));   
   RLP_Mag_Suc_Signal_dFF(k,:)=Signals_working_dFF(RLP_Mag_Suc(k,1),:);    
    end
end


%Aligning press-entry sequence signal data to mag entries 
if length(LLP_Mag_Pel)>0;   %execute operation only for LPel animals
for k=1:length(First_Mag_Pel);
    if First_Mag_Pel(k,1)<= length(Signals_working_dFF(:,1));  
   magaligned_LLP_Mag_Pel_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Pel(k,1),:);    
    end
end
end

if length(LLP_Mag_Suc)>0;   %execute operation only for LSuc animals
for k=1:length(First_Mag_Suc);
    if First_Mag_Suc(k,1)<= length(Signals_working_dFF(:,1));   
   magaligned_LLP_Mag_Suc_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Suc(k,1),:);    
    end
end
end


for k=1:length(magaligned_LLP_Mag_Empty);
    if magaligned_LLP_Mag_Empty(k,1)<= length(Signals_working_dFF(:,1)); 
   magaligned_LLP_Mag_Empty_Signal_dFF(k,:)=Signals_working_dFF(magaligned_LLP_Mag_Empty(k,1),:);    
    end
end

if length(RLP_Mag_Pel)>0;   %execute operation only for RPel animals
for k=1:length(First_Mag_Pel);
    if First_Mag_Pel(k,1)<= length(Signals_working_dFF(:,1));
   magaligned_RLP_Mag_Pel_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Pel(k,1),:);    
    end
end
end

if length(RLP_Mag_Suc)>0;   %execute operation only for RSuc animals
for k=1:length(First_Mag_Suc);
    if First_Mag_Suc(k,1)<= length(Signals_working_dFF(:,1));   
   magaligned_RLP_Mag_Suc_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Suc(k,1),:);    
    end
end
end


for k=1:length(magaligned_RLP_Mag_Empty);
    if magaligned_RLP_Mag_Empty(k,1)<= length(Signals_working_dFF(:,1)); 
   magaligned_RLP_Mag_Empty_Signal_dFF(k,:)=Signals_working_dFF(magaligned_RLP_Mag_Empty(k,1),:);    
    end
end



%storing signal data for all event sequence types in cell arrays
LLP_LLP_Signal_all{h,x}=LLP_LLP_Signal_dFF;
LLP_Mag_Empty_Signal_all{h,x}=LLP_Mag_Empty_Signal_dFF;
LLP_Mag_Pel_Signal_all{h,x}=LLP_Mag_Pel_Signal_dFF;
LLP_Mag_Suc_Signal_all{h,x}=LLP_Mag_Suc_Signal_dFF;
RLP_RLP_Signal_all{h,x}=RLP_RLP_Signal_dFF;
RLP_Mag_Empty_Signal_all{h,x}=RLP_Mag_Empty_Signal_dFF;
RLP_Mag_Pel_Signal_all{h,x}=RLP_Mag_Pel_Signal_dFF;
RLP_Mag_Suc_Signal_all{h,x}=RLP_Mag_Suc_Signal_dFF;
magaligned_LLP_Mag_Pel_Signal_all{h,x}=magaligned_LLP_Mag_Pel_Signal_dFF;
magaligned_LLP_Mag_Suc_Signal_all{h,x}=magaligned_LLP_Mag_Suc_Signal_dFF;
magaligned_LLP_Mag_Empty_Signal_all{h,x}=magaligned_LLP_Mag_Empty_Signal_dFF;
magaligned_RLP_Mag_Pel_Signal_all{h,x}=magaligned_RLP_Mag_Pel_Signal_dFF;
magaligned_RLP_Mag_Suc_Signal_all{h,x}=magaligned_RLP_Mag_Suc_Signal_dFF;
magaligned_RLP_Mag_Empty_Signal_all{h,x}=magaligned_RLP_Mag_Empty_Signal_dFF;

    end

end

end
