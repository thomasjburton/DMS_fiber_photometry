

%to categorize lever presses according to what preceded them- another
%press, mag entry or nothing (First press), and get signal data for each presstype.



LLP=1; %left lever press
RLP=4; %right lever press   
PEL=16; %pellet reward
SUC=8; %sucrose reward
Mag=2; %magazine entry


% Example given for cell array "PrtA_cells"


% PrtA
for id=1:length(PrtA_cells)
    h=PrtA_cells(id,1);

Events_timer=[];
Events_timer=Epocs_working_all{h,1};
Events_timer(:,[1,4])=[]; %create matrix with all events and their codes
Events_timer(:,4)=zeros; %Create columns for inserting relevant event timestamps
Events_timer(:,5)=zeros;
Events_timer(:,6)=zeros;
Events_timer(:,7)=zeros;


for i=1:length(Events_timer(:,3));  %remove negative event code rows (i.e. event terminations)
a=find (Events_timer(:,3)<1);
Events_timer(a,:)=[];
end

for i=2:length(Events_timer);  %insert mag entry events into col 4 and timestapms into col 5
Lia=[];
Lia=ismember(Events_timer(:,1),Locate_Mag_all{h,1});
Events_timer(:,8)=Lia;  
    if Events_timer(i,8)==1;
       Events_timer(i,4)=Mag;
       Events_timer(i,5)=Events_timer(i,2);
    else
       Events_timer(i,5)=Events_timer((i-1),5); 
    end
end

for i=2:length(Events_timer);  %insert LLP events into col 4 and timestapms into col 6
Lia=[];
Lia=ismember(Events_timer(:,1),Locate_LLP_all{h,1});
Events_timer(:,8)=Lia;  
    if Events_timer(i,8)==1;
       Events_timer(i,4)=LLP;
       Events_timer(i,6)=Events_timer(i,2);
    else
       Events_timer(i,6)=Events_timer((i-1),6); 
    end
end

for i=2:length(Events_timer);  %insert RLP events into col 4 and timestapms into col 7
Lia=[];
Lia=ismember(Events_timer(:,1),Locate_RLP_all{h,1});
Events_timer(:,8)=Lia;  
    if Events_timer(i,8)==1;
       Events_timer(i,4)=RLP;
       Events_timer(i,7)=Events_timer(i,2);
    else
       Events_timer(i,7)=Events_timer((i-1),7); 
    end
end

% Find ITI, SD & Mean for LLP & RLP

ITI_LLP=[];
ITI_RLP=[];

for a=2:size (Events_timer) 
ITI_LLP(a,1)=Events_timer(a,6)-Events_timer(a-1,6);
end
ITI_LLP=nonzeros(ITI_LLP);

ITI_LLP(1,1)=0;
b=find(ITI_LLP>300);  %exclude inter-press-intervals >5 minutes
ITI_LLP(b)=0;

ITI_LLP=nonzeros(ITI_LLP);

ITI_LLP_mean=mean(ITI_LLP);
ITI_LLP_std=std(ITI_LLP);

ITI_LLP_std_all{h,1}=ITI_LLP_std;
ITI_LLP_mean_all{h,1}=ITI_LLP_mean;

% Find ITI, SD & Mean for RLP
for a=2:size (Events_timer) 
ITI_RLP(a,1)=Events_timer(a,7)-Events_timer(a-1,7);
end
ITI_RLP=nonzeros(ITI_RLP);

ITI_RLP(1,1)=0;
c=find(ITI_RLP>300);
ITI_RLP(c)=0;

ITI_RLP=nonzeros(ITI_RLP);

ITI_RLP_mean=mean(ITI_RLP);
ITI_RLP_std=std(ITI_RLP);

ITI_RLP_std_all{h,1}=ITI_RLP_std;
ITI_RLP_mean_all{h,1}=ITI_RLP_mean;

Events_timer(:,8)=[];
LLP_pre_LLP_nobase=[];
Mag_pre_LLP_nobase=[];
Nil_pre_LLP_nobase=[];

for i=2:length(Events_timer); %categorizing LLP according to preceding event
    if Events_timer(i,4)==LLP & (Events_timer(i,6)-Events_timer(i,5))>(Events_timer(i,6)-Events_timer(i-1,6))& (Events_timer(i,6)-Events_timer(i-1,6))<ITI_LLP_std;
    LLP_pre_LLP_nobase(i,1)=Events_timer(i,1); 
    elseif Events_timer(i,4)==LLP & (Events_timer(i,6)-Events_timer(i,5))<(Events_timer(i,6)-Events_timer((i-1),6))& (Events_timer(i,6)-Events_timer(i,5))<ITI_LLP_std;
    Mag_pre_LLP_nobase(i,1)=Events_timer(i,1); 
    elseif Events_timer(i,4)==LLP & (Events_timer(i,6)-Events_timer(i,5))>ITI_LLP_std & (Events_timer(i,6)-Events_timer((i-1),6))>ITI_LLP_std;
    Nil_pre_LLP_nobase(i,1)=Events_timer(i,1);
    end
end

RLP_pre_RLP_nobase=[];
Mag_pre_RLP_nobase=[];
Nil_pre_RLP_nobase=[];

for i=2:length(Events_timer); %categorizing RLP according to preceding event
    if Events_timer(i,4)==RLP & (Events_timer(i,7)-Events_timer(i,5))>(Events_timer(i,7)-Events_timer(i-1,7))& (Events_timer(i,7)-Events_timer(i-1,7))<ITI_RLP_std;
    RLP_pre_RLP_nobase(i,1)=Events_timer(i,1); 
    elseif Events_timer(i,4)==RLP & (Events_timer(i,7)-Events_timer(i,5))<(Events_timer(i,7)-Events_timer((i-1),7))& (Events_timer(i,7)-Events_timer(i,5))<ITI_RLP_std;
    Mag_pre_RLP_nobase(i,1)=Events_timer(i,1); 
    elseif Events_timer(i,4)==RLP & (Events_timer(i,7)-Events_timer(i,5))>ITI_RLP_std & (Events_timer(i,7)-Events_timer((i-1),7))>ITI_RLP_std;
    Nil_pre_RLP_nobase(i,1)=Events_timer(i,1);
    end
end

RLP_pre_RLP_nobase=nonzeros(RLP_pre_RLP_nobase);
Mag_pre_RLP_nobase=nonzeros(Mag_pre_RLP_nobase);
Nil_pre_RLP_nobase=nonzeros(Nil_pre_RLP_nobase);

LLP_pre_LLP_nobase=nonzeros(LLP_pre_LLP_nobase);
Mag_pre_LLP_nobase=nonzeros(Mag_pre_LLP_nobase);
Nil_pre_LLP_nobase=nonzeros(Nil_pre_LLP_nobase);

Nil_pre_RLP_nobase_all{h,1}=Nil_pre_RLP_nobase;
Mag_pre_RLP_nobase_all{h,1}=Mag_pre_RLP_nobase;
RLP_pre_RLP_nobase_all{h,1}=RLP_pre_RLP_nobase;
Nil_pre_LLP_nobase_all{h,1}=Nil_pre_LLP_nobase;
Mag_pre_LLP_nobase_all{h,1}=Mag_pre_LLP_nobase;
LLP_pre_LLP_nobase_all{h,1}=LLP_pre_LLP_nobase;
Events_timer_all{h,1}=Events_timer;

%excluding presses where preceding events fall within the baseline

LLP_base_LLP=[];

for i=2:length(Events_timer); %categorizing LLP according to preceding event
    if Events_timer(i,4)==LLP & (Events_timer(i,6)-Events_timer(i-1,6))<5 & (Events_timer(i,6)-Events_timer(i-1,6))>2.5 | Events_timer(i,4)==LLP & (Events_timer(i,6)-Events_timer(i,5))<5 & (Events_timer(i,6)-Events_timer(i,5))>2.5;
        LLP_base_LLP(i,1)=Events_timer(i,1);   
    end
end

RLP_base_RLP=[];

for i=2:length(Events_timer); %categorizing RLP according to preceding event
    if Events_timer(i,4)==RLP & (Events_timer(i,7)-Events_timer(i-1,7))<5 & (Events_timer(i,7)-Events_timer(i-1,7))>2.5|Events_timer(i,4)==RLP & (Events_timer(i,7)-Events_timer(i,5))<5 & (Events_timer(i,7)-Events_timer(i,5))>2.5; ;
    RLP_base_RLP(i,1)=Events_timer(i,1); 
    end
end

LLP_base_LLP=nonzeros(LLP_base_LLP);

RLP_base_RLP=nonzeros(RLP_base_RLP);

LLP_pre_LLP=setdiff(LLP_pre_LLP_nobase,LLP_base_LLP);
Mag_pre_LLP=setdiff(Mag_pre_LLP_nobase,LLP_base_LLP);
Nil_pre_LLP=setdiff(Nil_pre_LLP_nobase,LLP_base_LLP);
RLP_pre_RLP=setdiff(RLP_pre_RLP_nobase,RLP_base_RLP);
Mag_pre_RLP=setdiff(Mag_pre_RLP_nobase,RLP_base_RLP);
Nil_pre_RLP=setdiff(Nil_pre_RLP_nobase,RLP_base_RLP);

LLP_pre_LLP_all{h,1}=LLP_pre_LLP;
Mag_pre_LLP_all{h,1}=Mag_pre_LLP;
Nil_pre_LLP_all{h,1}=Nil_pre_LLP;
RLP_pre_RLP_all{h,1}=RLP_pre_RLP;
Mag_pre_RLP_all{h,1}=Mag_pre_RLP;
Nil_pre_RLP_all{h,1}=Nil_pre_RLP;

%Extracting dFF signals
Signals_working_dFF=[];
Signals_working_dFF=local_2s_z_dFF_std_all{h,1}; 

%extraction of signal data under broad categories of behavioural data
LLP_pre_LLP_Signal_dFF=[];
Mag_pre_LLP_Signal_dFF=[];
Nil_pre_LLP_Signal_dFF=[];

for k=1:length(LLP_pre_LLP);
    if LLP_pre_LLP(k,1)< length(Signals_working_dFF(:,1)); 
   LLP_pre_LLP_Signal_dFF(k,:)=Signals_working_dFF(LLP_pre_LLP(k,1),:);    
    end
end

for k=1:length(Mag_pre_LLP);
   if Mag_pre_LLP(k,1)< length(Signals_working_dFF(:,1)); 
     Mag_pre_LLP_Signal_dFF(k,:)=Signals_working_dFF(Mag_pre_LLP(k,1),:);    
   end
end
    
for k=1:length(Nil_pre_LLP);
   if Nil_pre_LLP(k,1)< length(Signals_working_dFF(:,1)); 
     Nil_pre_LLP_Signal_dFF(k,:)=Signals_working_dFF(Nil_pre_LLP(k,1),:);    
   end
end

LLP_pre_LLP_Signal_all{h,1}=LLP_pre_LLP_Signal_dFF;
Mag_pre_LLP_Signal_all{h,1}=Mag_pre_LLP_Signal_dFF;
Nil_pre_LLP_Signal_all{h,1}=Nil_pre_LLP_Signal_dFF;

%extraction of signal data under broad categories of behavioural data
RLP_pre_RLP_Signal_dFF=[];
Mag_pre_RLP_Signal_dFF=[];
Nil_pre_RLP_Signal_dFF=[];

for k=1:length(RLP_pre_RLP);
    if RLP_pre_RLP(k,1)< length(Signals_working_dFF(:,1)); 
   RLP_pre_RLP_Signal_dFF(k,:)=Signals_working_dFF(RLP_pre_RLP(k,1),:);    
    end
end

for k=1:length(Mag_pre_RLP);
   if Mag_pre_RLP(k,1)< length(Signals_working_dFF(:,1)); 
     Mag_pre_RLP_Signal_dFF(k,:)=Signals_working_dFF(Mag_pre_RLP(k,1),:);    
   end
end
    
for k=1:length(Nil_pre_RLP)
   if Nil_pre_RLP(k,1)< length(Signals_working_dFF(:,1)); 
     Nil_pre_RLP_Signal_dFF(k,:)=Signals_working_dFF(Nil_pre_RLP(k,1),:);    
   end
end

RLP_pre_RLP_Signal_all{h,1}=RLP_pre_RLP_Signal_dFF;
Mag_pre_RLP_Signal_all{h,1}=Mag_pre_RLP_Signal_dFF;
Nil_pre_RLP_Signal_all{h,1}=Nil_pre_RLP_Signal_dFF;

end


%//////////// AUC calculations///////////

%Given here for one example dataset "CONTRA_dFF_all"

bins=length(CONTRA_dFF_all(1,:));   %choose the variable that you want to calculate AUC for. Note, you will need to use cell2mat if organised in cellarray form
bps=bins/10;    %bins per s
epoc_on=bins/2;

EOI(1,1)=round(epoc_on-0.2*bps);       
EOI(1,2)=round(epoc_on+0.2*bps);

AUC_CONTRA=[];

%AUC 
for i=1:length(CONTRA_dFF_all(:,1));
AUC_CONTRA(i,1)=trapz(CONTRA_dFF_all(i,EOI(1,1):EOI(1,2)));
end;


