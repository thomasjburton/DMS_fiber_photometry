% To be used after main event extraction, uses pre-extracted events to
% simulate reward learning using a basic error rule
%Example below for rats with left lever-pellet (right lever-sucrose),
%create cell array "Lpel" with each cell containing an array for each animal, sessions in order.

Pel_Rew=18;
Suc_Rew=10;
LLP_No_Rew=3;
RLP_No_Rew=6;
LLP=1;
RLP=4;
alpha=0.3;
lambdapos=100;
lambdaneg=0;
startVpel=0;
startVsuc=0;


%PrtA L-pel R-suc
for i = 1:length(Lpel);
    b=[];
     b=Lpel{i};
for id=1:length(b)
    h=b(1,id);

RW_Epocs=[];
RW_Epocs=Epocs_working_all{h,1};
RW_Epocs(:,[1,3,4])=[]; %create matrix with all events and their codes
RW_Epocs(:,3)=zeros; %Create 3rd column for inserting relevant event combination codes

    
for i=1:length(RW_Epocs(:,2));  %remove negative event code rows (i.e. event terminations)
a=find (RW_Epocs(:,2)<1);
RW_Epocs(a,:)=[];
end


for i=1:length(RW_Epocs);  %insert (earned) pel retrieval events into col 4
Lia=[];
Lia=ismember(RW_Epocs(:,1),First_Mag_Pel_all{h,1});
RW_Epocs(:,4)=Lia;  
    if RW_Epocs(i,4)==1;
       RW_Epocs(i,3)=Pel_Rew;
    end
end

for i=1:length(RW_Epocs);  %insert (earned) suc retrieval events into col 4
Lia=[];
Lia=ismember(RW_Epocs(:,1),First_Mag_Suc_all{h,1});
RW_Epocs(:,4)=Lia;  
    if RW_Epocs(i,4)==1;
       RW_Epocs(i,3)=Suc_Rew;
    end
end

for i=1:length(RW_Epocs);  %insert empty mags after LLP events into col 4
Lia=[];
Lia=ismember(RW_Epocs(:,1),magaligned_LLP_Mag_Empty_all{h,1});
RW_Epocs(:,4)=Lia;  
    if RW_Epocs(i,4)==1;
       RW_Epocs(i,3)=LLP_No_Rew;
    end
end

for i=1:length(RW_Epocs);  %insert empty mags after RLP events into col 4
Lia=[];
Lia=ismember(RW_Epocs(:,1),magaligned_RLP_Mag_Empty_all{h,1});
RW_Epocs(:,4)=Lia;  
    if RW_Epocs(i,4)==1;
       RW_Epocs(i,3)=RLP_No_Rew;
    end
end


lambdacolpel=[];%building columns for lambda values for each event type
lambdacolsuc=[];

lambdacolpel=RW_Epocs(:,3);
lambdacolpel(lambdacolpel==0)=NaN;
lambdacolpel(lambdacolpel==Suc_Rew)=NaN;
lambdacolpel(lambdacolpel==RLP_No_Rew)=NaN;
lambdacolpel(lambdacolpel==Pel_Rew)=lambdapos;
lambdacolpel(lambdacolpel==LLP_No_Rew)=lambdaneg;

lambdacolsuc=RW_Epocs(:,3);
lambdacolsuc(lambdacolsuc==0)=NaN;
lambdacolsuc(lambdacolsuc==Pel_Rew)=NaN;
lambdacolsuc(lambdacolsuc==LLP_No_Rew)=NaN;
lambdacolsuc(lambdacolsuc==Suc_Rew)=lambdapos;
lambdacolsuc(lambdacolsuc==RLP_No_Rew)=lambdaneg;

RW_Epocs(:,4)=lambdacolpel;
RW_Epocs(:,5)=lambdacolsuc;

RW_Epocs(:,6)=zeros; %building column for deltaV for each pel event
RW_Epocs(:,7)=zeros;    %column for sumV for pel
RW_Epocs(1,7)=startVpel;

RW_Epocs(:,8)=zeros; %building column for deltaV for each suc event
RW_Epocs(:,9)=zeros;    %column for sumV for suc
RW_Epocs(1,9)=startVsuc;

%RW calculation for L-Pel

for i=2:length(RW_Epocs);
    if RW_Epocs(i,4)==lambdapos;
 RW_Epocs(i,6)=alpha.*((RW_Epocs(i,4))-(RW_Epocs(i-1,7)));
    elseif RW_Epocs(i,4)==lambdaneg;
 RW_Epocs(i,6)=alpha.*((RW_Epocs(i,4))-(RW_Epocs(i-1,7)));   
    end
    RW_Epocs(i,7)=(RW_Epocs(i-1,7)) + (RW_Epocs(i,6));
end

%RW calculation for R-Suc

for i=2:length(RW_Epocs);
    if RW_Epocs(i,5)==lambdapos;
 RW_Epocs(i,8)=alpha.*((RW_Epocs(i,5))-(RW_Epocs(i-1,9)));
    elseif RW_Epocs(i,5)==lambdaneg;
 RW_Epocs(i,8)=alpha.*((RW_Epocs(i,5))-(RW_Epocs(i-1,9)));   
    end
    RW_Epocs(i,9)=(RW_Epocs(i-1,9)) + (RW_Epocs(i,8));
end
 
startVpel=RW_Epocs(end,7);
startVsuc=RW_Epocs(end,9);

RW_Epocs_all{h,1}=RW_Epocs;
end

startVpel=0;
startVsuc=0;

end