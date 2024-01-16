
RLP=4;
PEL=16;
SUC=8;
ENTRY_CODE=2;

for id=1:length(PrtA_cells)
    h=PrtA_cells(id,1);
    Epocs_working=[];
    Epocs_working=EpocCat_allanimals{h,1};
    Epocs_working(:,5)=zeros;          
    
    
    %creating the on/off record by subtracting epoc(n-1) from epoch(n)
     for i=2:length(Epocs_working);
        Epocs_working(i,5)=Epocs_working(i,1)-Epocs_working(i-1,1);
     end

    for i=1:length(Epocs_working);
        Epocs_working(i,2)=i;
    end
    
    Epocs_working_all{h,1}=Epocs_working;
    
   
    %locating broad/elemetal events of interest. Must be searching the on/off record calculated above       
    Locate_LLP=[];
    Locate_RLP=[];
    Locate_Mag=[];
    Locate_Pel=[];
    Locate_Suc=[];
    Locate_Mag=[];
    Locate_LLP=find(Epocs_working(:,5)==LLP | Epocs_working(:,5)==LLP+PEL | Epocs_working(:,5)==LLP+SUC) ; 
    Locate_RLP=find(Epocs_working(:,5)==RLP| Epocs_working(:,5)==RLP+PEL | Epocs_working(:,5)==RLP+SUC) ;
    Locate_Mag=find(Epocs_working(:,5)==ENTRY_CODE) ;
    Locate_Pel=find(Epocs_working(:,5)==PEL | Epocs_working(:,5)==LLP+PEL | Epocs_working(:,5)==RLP+PEL) ;
    Locate_Suc=find(Epocs_working(:,5)==SUC | Epocs_working(:,5)==LLP+SUC | Epocs_working(:,5)==RLP+SUC) ;

    Mag_working=[];
    for i=1:length(Epocs_working);
    if Epocs_working(i,5) == 2;         
        Mag_working(i,1:5)=Epocs_working(i,1:5);
    else Mag_working(i,1:5)=0;
    end 
    end
    
    Locate_LLP_all{h,1}=Locate_LLP;
    Locate_RLP_all{h,1}=Locate_RLP;
    Locate_Mag_all{h,1}=Locate_Mag;
    
 
    %categorising outcomes as earned or free    
Earned_O=[];
Free_O=[];
Locate_Pel_Earned=[];
Locate_Pel_Free=[];
Locate_Pel_Earned=[];
Locate_Pel_Free=[];
Locate_Earned=[];
Locate_Free=[];

   for i=2:length(Epocs_working);
   if Epocs_working(i,5)==PEL & Epocs_working(i-1,5)==LLP | Epocs_working(i,5)==PEL & Epocs_working(i-1,5)==RLP | Epocs_working(i,5)==SUC & Epocs_working(i-1,5)==LLP | Epocs_working(i,5)==SUC & Epocs_working(i-1,5)==RLP | Epocs_working(i,5)==PEL+LLP | Epocs_working(i,5)==PEL+RLP | Epocs_working(i,5)==SUC+LLP | Epocs_working(i,5)==SUC+RLP | Epocs_working(i,5)==PEL &(find(Epocs_working(:,5)==-LLP & Epocs_working(:,2)>i,1))<Locate_LLP(find(Locate_LLP(:,1)>i,1)) | Epocs_working(i,5)==PEL &(find(Epocs_working(:,5)==-RLP & Epocs_working(:,2)>i,1))<Locate_RLP(find(Locate_RLP(:,1)>i,1)) | Epocs_working(i,5)==SUC &(find(Epocs_working(:,5)==-RLP & Epocs_working(:,2)>i,1))<Locate_RLP(find(Locate_RLP(:,1)>i,1)) | Epocs_working(i,5)==SUC &(find(Epocs_working(:,5)==-RLP & Epocs_working(:,2)>i,1))<Locate_RLP(find(Locate_RLP(:,1)>i,1))
   
       Earned_O(i,1:5)=Epocs_working(i,1:5);
   elseif Epocs_working(i,5)==PEL | Epocs_working(i,5)==SUC
       Free_O(i,1:5)=Epocs_working(i,1:5);
               end
   end

   Locate_Pel_Earned=find(Earned_O(:,5)==PEL | Earned_O(:,5)==PEL+LLP | Earned_O(:,5)==PEL+RLP);    %changed all (:,1) to (:,5)
   Locate_Pel_Free=find(Free_O(:,5)==PEL | Free_O(:,5)==PEL+LLP | Free_O(:,5)==PEL+RLP | Free_O(:,5)==PEL+ENTRY_CODE); %added in PEL+ENTRY_CODE
   Locate_Suc_Earned=find(Earned_O(:,5)==SUC | Earned_O(:,5)==SUC+LLP | Earned_O(:,5)==SUC+RLP);
   Locate_Suc_Free=find(Free_O(:,5)==SUC | Free_O(:,5)==SUC+LLP | Free_O(:,5)==SUC+RLP | Free_O(:,5)==SUC+ENTRY_CODE);
   Locate_Earned=nonzeros(Earned_O(:,2));
   Locate_Free=nonzeros(Free_O(:,2));
   
%Locating Mag entries that putatively contain outcomes
  First_Mag_Pel=[];
  First_Mag_Suc=[];
  First_Mag_Pel_Earned=[];
  First_Mag_Pel_Free=[];
  First_Mag_Suc_Earned=[];
  First_Mag_Suc_Free=[];
  First_Mag_Earned=[];
  First_Mag_Free=[];
 
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
  
   for j=1:size(Locate_Pel_Earned)-1;
       if Locate_Pel_Earned(j)<Locate_Mag(end) & Locate_Pel_Earned(j+1)>find(Mag_working(:,2)>Locate_Pel_Earned(j),1);
    First_Mag_Pel_Earned(j,1)=find(Mag_working(:,2)>Locate_Pel_Earned(j),1);
       end 
   end
   if length(Locate_Pel_Earned)>0;
   j=length(Locate_Pel_Earned);
   if Locate_Pel_Earned(j)<Locate_Mag(end)
    First_Mag_Pel_Earned(j,1)=find(Mag_working(:,2)>Locate_Pel_Earned(j),1);
   end
  First_Mag_Pel_Earned=nonzeros(First_Mag_Pel_Earned);
   end
   
   for j=1:size(Locate_Pel_Free)-1;
       if Locate_Pel_Free(j)<Locate_Mag(end) & Locate_Pel_Free(j+1)>find(Mag_working(:,2)>Locate_Pel_Free(j),1);
    First_Mag_Pel_Free(j,1)=find(Mag_working(:,2)>Locate_Pel_Free(j),1);
       end 
   end
   if length(Locate_Pel_Free)>0;
   j=length(Locate_Pel_Free);
   if Locate_Pel_Free(j)<Locate_Mag(end)
    First_Mag_Pel_Free(j,1)=find(Mag_working(:,2)>Locate_Pel_Free(j),1);
   end
  First_Mag_Pel_Free=nonzeros(First_Mag_Pel_Free);
   end
   
   for j=1:size(Locate_Suc_Earned)-1;
       if Locate_Suc_Earned(j)<Locate_Mag(end) & Locate_Suc_Earned(j+1)>find(Mag_working(:,2)>Locate_Suc_Earned(j),1);
    First_Mag_Suc_Earned(j,1)=find(Mag_working(:,2)>Locate_Suc_Earned(j),1);
       end 
   end
   if length(Locate_Suc_Earned)>0;
   j=length(Locate_Suc_Earned);
   if Locate_Suc_Earned(j)<Locate_Mag(end)
    First_Mag_Suc_Earned(j,1)=find(Mag_working(:,2)>Locate_Suc_Earned(j),1);
   end
  First_Mag_Suc_Earned=nonzeros(First_Mag_Suc_Earned);
   end
   
   for j=1:size(Locate_Suc_Free)-1;
       if Locate_Suc_Free(j)<Locate_Mag(end) & Locate_Suc_Free(j+1)>find(Mag_working(:,2)>Locate_Suc_Free(j),1);
    First_Mag_Suc_Free(j,1)=find(Mag_working(:,2)>Locate_Suc_Free(j),1);
       end 
   end
   if length(Locate_Suc_Free)>0;
   j=length(Locate_Suc_Free);
   if Locate_Suc_Free(j)<Locate_Mag(end)
    First_Mag_Suc_Free(j,1)=find(Mag_working(:,2)>Locate_Suc_Free(j),1);
   end
  First_Mag_Suc_Free=nonzeros(First_Mag_Suc_Free);
   end
   
   for j=1:size(Locate_Earned)-1;
       if Locate_Earned(j)<Locate_Mag(end) & Locate_Earned(j+1)>find(Mag_working(:,2)>Locate_Earned(j),1);
    First_Mag_Earned(j,1)=find(Mag_working(:,2)>Locate_Earned(j),1);
       end 
   end
   if length(Locate_Earned)>0;
   j=length(Locate_Earned);
   if Locate_Earned(j)<Locate_Mag(end)    
    First_Mag_Earned(j,1)=find(Mag_working(:,2)>Locate_Earned(j),1);       
   end
  First_Mag_Earned=nonzeros(First_Mag_Earned);
   end
   
   for j=1:size(Locate_Free)-1;
       if Locate_Free(j)<Locate_Mag(end) & Locate_Free(j+1)>find(Mag_working(:,2)>Locate_Free(j),1);
    First_Mag_Free(j,1)=find(Mag_working(:,2)>Locate_Free(j),1);
       end 
   end
   if length(Locate_Free)>0;
   j=length(Locate_Free);
   if Locate_Free(j)<Locate_Mag(end)   
    First_Mag_Free(j,1)=find(Mag_working(:,2)>Locate_Free(j),1);       
   end
  First_Mag_Free=nonzeros(First_Mag_Free);
   end
   
 First_Mag_Pel_all{h,1}=First_Mag_Pel;
 First_Mag_Suc_all{h,1}=First_Mag_Suc;
 First_Mag_Pel_Earned_all{h,1}=First_Mag_Pel_Earned;
 First_Mag_Suc_Earned_all{h,1}=First_Mag_Suc_Earned;
 First_Mag_Pel_Free_all{h,1}=First_Mag_Pel_Free;
 First_Mag_Suc_Free_all{h,1}=First_Mag_Suc_Free;
 First_Mag_Earned_all{h,1}=First_Mag_Earned;
 First_Mag_Free_all{h,1}=First_Mag_Free;

%locating empty mag entries
Locate_Mag_Empty=[];
Locate_Mag_Empty=setdiff(Locate_Mag,First_Mag_Pel);
Locate_Mag_Empty=setdiff(Locate_Mag_Empty,First_Mag_Suc);

Locate_Mag_Empty_all{h,1}=Locate_Mag_Empty;

%Extracting dFF signals
Signals_working_dFF=[];
Signals_working_dFF=local_2s_z_dFF_std_all{h,1}; 

First_Mag_Pel_Signal_dFF=[];
First_Mag_Suc_Signal_dFF=[];
Mag_Empty_Signal_dFF=[];

First_Mag_Pel_Earned_Signal_dFF=[];
First_Mag_Pel_Free_Signal_dFF=[];
First_Mag_Suc_Earned_Signal_dFF=[];
First_Mag_Suc_Free_Signal_dFF=[];


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
    
for k=1:length(Locate_Mag_Empty);
   if Locate_Mag_Empty(k,1)<= length(Signals_working_dFF(:,1)); 
   Mag_Empty_Signal_dFF(k,:)=Signals_working_dFF(Locate_Mag_Empty(k,1),:);    
   end
end

for k=1:length(First_Mag_Pel_Earned);
    if First_Mag_Pel_Earned(k,1)<= length(Signals_working_dFF(:,1));
   First_Mag_Pel_Earned_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Pel_Earned(k,1),:);       
    end
end

for k=1:length(First_Mag_Pel_Free);
    if First_Mag_Pel_Free(k,1)<= length(Signals_working_dFF(:,1));
   First_Mag_Pel_Free_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Pel_Free(k,1),:);    
    end
end

for k=1:length(First_Mag_Suc_Earned);
    if First_Mag_Suc_Earned(k,1)<= length(Signals_working_dFF(:,1));
   First_Mag_Suc_Earned_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Suc_Earned(k,1),:);        
    end
end

for k=1:length(First_Mag_Suc_Free);
    if First_Mag_Suc_Free(k,1)<= length(Signals_working_dFF(:,1));
   First_Mag_Suc_Free_Signal_dFF(k,:)=Signals_working_dFF(First_Mag_Suc_Free(k,1),:);    
    end
end


First_Mag_Pel_Signal_all{h,1}=First_Mag_Pel_Signal_dFF;

First_Mag_Suc_Signal_all{h,1}=First_Mag_Suc_Signal_dFF;

Mag_Empty_Signal_all{h,1}=Mag_Empty_Signal_dFF;

First_Mag_Pel_Earned_Signal_all{h,1}=First_Mag_Pel_Earned_Signal_dFF;

First_Mag_Pel_Free_Signal_all{h,1}=First_Mag_Pel_Free_Signal_dFF;

First_Mag_Suc_Earned_Signal_all{h,1}=First_Mag_Suc_Earned_Signal_dFF;

First_Mag_Suc_Free_Signal_all{h,1}=First_Mag_Suc_Free_Signal_dFF;



%categorising LPs
LLP_LLP=[];
LLP_Mag=[];
LLP_Mag_Empty=[];
LLP_Mag_Pel_Earned=[];
LLP_Mag_Suc_Earned=[];
LLP_Mag_Pel_Free=[];
LLP_Mag_Suc_Free=[];
RLP_RLP=[];
RLP_Mag=[];
RLP_Mag_Empty=[];
RLP_Mag_Pel_Earned=[];
RLP_Mag_Suc_Earned=[];
RLP_Mag_Pel_Free=[];
RLP_Mag_Suc_Free=[];


for j=1:size(Locate_LLP)-1; %not looking at final LP of session
       if Locate_LLP(j)<Locate_Mag(end) & Locate_LLP(j+1)<find(Mag_working(:,2)>Locate_LLP(j),1);
    LLP_LLP(j,1)=Locate_LLP(j);
       end 
end

LLP_LLP=nonzeros(LLP_LLP);
LLP_Mag=setdiff(Locate_LLP,LLP_LLP);

LLP_LLP_all{h,1}=LLP_LLP;
LLP_Mag_all{h,1}=LLP_Mag;

for j=1:size(LLP_Mag); %not looking at final LP of session
    M_O_working=zeros(5,2); %create a working variable that will locate the first occurence of a mag entry type for a give LP-Mag combination
    
     if length(Locate_Mag_Empty)>0;
    if LLP_Mag(j)<Locate_Mag_Empty(end);
    M_O_working(1,1)=find(Locate_Mag_Empty(:,1)>LLP_Mag(j),1);
    M_O_working(1,2)=Locate_Mag_Empty(M_O_working(1,1),1);
    end
    end
    
    if length(First_Mag_Suc_Earned)>0;
    if LLP_Mag(j)<First_Mag_Suc_Earned(end);
    M_O_working(2,1)=find(First_Mag_Suc_Earned(:,1)>LLP_Mag(j),1);
    M_O_working(2,2)=First_Mag_Suc_Earned(M_O_working(2,1),1);
    end
    end
   
    if length(First_Mag_Pel_Earned)>0;
    if LLP_Mag(j)<First_Mag_Pel_Earned(end);
    M_O_working(3,1)=find(First_Mag_Pel_Earned(:,1)>LLP_Mag(j),1);
    M_O_working(3,2)=First_Mag_Pel_Earned(M_O_working(3,1),1);
    end
    end
    
    if length(First_Mag_Suc_Free)>0;
    if LLP_Mag(j)<First_Mag_Suc_Free(end);
    M_O_working(4,1)=find(First_Mag_Suc_Free(:,1)>LLP_Mag(j),1);
    M_O_working(4,2)=First_Mag_Suc_Free(M_O_working(4,1),1);
    end
    end
    
    if length(First_Mag_Pel_Free)>0;
    if LLP_Mag(j)<First_Mag_Pel_Free(end);
    M_O_working(5,1)=find(First_Mag_Pel_Free(:,1)>LLP_Mag(j),1);
    M_O_working(5,2)=First_Mag_Pel_Free(M_O_working(5,1),1);
    end
    end
    
    if M_O_working(1,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(1,2)>0;
    LLP_Mag_Empty(j,1)=LLP_Mag(j,1);
    elseif M_O_working(2,2)== min(nonzeros(M_O_working(:,2))) & M_O_working(2,2)>0;
    LLP_Mag_Suc_Earned(j,1)=LLP_Mag(j,1);
    elseif M_O_working(3,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(3,2)>0;
    LLP_Mag_Pel_Earned(j,1)=LLP_Mag(j,1);
    elseif M_O_working(4,2)== min(nonzeros(M_O_working(:,2))) & M_O_working(4,2)>0;
    LLP_Mag_Suc_Free(j,1)=LLP_Mag(j,1);
    elseif M_O_working(5,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(5,2)>0;
    LLP_Mag_Pel_Free(j,1)=LLP_Mag(j,1);
    end
end

LLP_Mag_Empty=nonzeros(LLP_Mag_Empty);
LLP_Mag_Pel_Earned=nonzeros(LLP_Mag_Pel_Earned);
LLP_Mag_Suc_Earned=nonzeros(LLP_Mag_Suc_Earned);
LLP_Mag_Pel_Free=nonzeros(LLP_Mag_Pel_Free);
LLP_Mag_Suc_Free=nonzeros(LLP_Mag_Suc_Free);

LLP_Mag_Empty_all{h,1}=LLP_Mag_Empty;
LLP_Mag_Pel_Earned_all{h,1}=LLP_Mag_Pel_Earned;
LLP_Mag_Suc_Earned_all{h,1}=LLP_Mag_Suc_Earned;
LLP_Mag_Pel_Free_all{h,1}=LLP_Mag_Pel_Free;
LLP_Mag_Suc_Free_all{h,1}=LLP_Mag_Suc_Free;


for j=1:size(Locate_RLP)-1; %not looking at final LP of session
       if Locate_RLP(j)<Locate_Mag(end) & Locate_RLP(j+1)<find(Mag_working(:,2)>Locate_RLP(j),1);
    RLP_RLP(j,1)=Locate_RLP(j);
       end 
end

RLP_RLP=nonzeros(RLP_RLP);
RLP_Mag=setdiff(Locate_RLP,RLP_RLP);

RLP_RLP_all{h,1}=RLP_RLP;
RLP_Mag_all{h,1}=RLP_Mag;

for j=1:size(RLP_Mag); %not looking at final LP of session
    M_O_working=zeros(5,2); %create a working variable that will locate the first occurence of a mag entry type for a give LP-Mag combination
    
     if length(Locate_Mag_Empty)>0;
    if RLP_Mag(j)<Locate_Mag_Empty(end);
    M_O_working(1,1)=find(Locate_Mag_Empty(:,1)>RLP_Mag(j),1);
    M_O_working(1,2)=Locate_Mag_Empty(M_O_working(1,1),1);
    end
    end
    
    if length(First_Mag_Suc_Earned)>0;
    if RLP_Mag(j)<First_Mag_Suc_Earned(end);
    M_O_working(2,1)=find(First_Mag_Suc_Earned(:,1)>RLP_Mag(j),1);
    M_O_working(2,2)=First_Mag_Suc_Earned(M_O_working(2,1),1);
    end
    end
   
    if length(First_Mag_Pel_Earned)>0;
    if RLP_Mag(j)<First_Mag_Pel_Earned(end);
    M_O_working(3,1)=find(First_Mag_Pel_Earned(:,1)>RLP_Mag(j),1);
    M_O_working(3,2)=First_Mag_Pel_Earned(M_O_working(3,1),1);
    end
    end
    
    if length(First_Mag_Suc_Free)>0;
    if RLP_Mag(j)<First_Mag_Suc_Free(end);
    M_O_working(4,1)=find(First_Mag_Suc_Free(:,1)>RLP_Mag(j),1);
    M_O_working(4,2)=First_Mag_Suc_Free(M_O_working(4,1),1);
    end
    end
    
    if length(First_Mag_Pel_Free)>0;
    if RLP_Mag(j)<First_Mag_Pel_Free(end);
    M_O_working(5,1)=find(First_Mag_Pel_Free(:,1)>RLP_Mag(j),1);
    M_O_working(5,2)=First_Mag_Pel_Free(M_O_working(5,1),1);
    end
    end
    
    if M_O_working(1,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(1,2)>0;
    RLP_Mag_Empty(j,1)=RLP_Mag(j,1);
    elseif M_O_working(2,2)== min(nonzeros(M_O_working(:,2))) & M_O_working(2,2)>0;
    RLP_Mag_Suc_Earned(j,1)=RLP_Mag(j,1);
    elseif M_O_working(3,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(3,2)>0;
    RLP_Mag_Pel_Earned(j,1)=RLP_Mag(j,1);
    elseif M_O_working(4,2)== min(nonzeros(M_O_working(:,2))) & M_O_working(4,2)>0;
    RLP_Mag_Suc_Free(j,1)=RLP_Mag(j,1);
    elseif M_O_working(5,2) == min(nonzeros(M_O_working(:,2))) & M_O_working(5,2)>0;
    RLP_Mag_Pel_Free(j,1)=RLP_Mag(j,1);
    end
end

RLP_Mag_Empty=nonzeros(RLP_Mag_Empty);
RLP_Mag_Pel_Earned=nonzeros(RLP_Mag_Pel_Earned);
RLP_Mag_Suc_Earned=nonzeros(RLP_Mag_Suc_Earned);
RLP_Mag_Pel_Free=nonzeros(RLP_Mag_Pel_Free);
RLP_Mag_Suc_Free=nonzeros(RLP_Mag_Suc_Free);

RLP_Mag_Empty_all{h,1}=RLP_Mag_Empty;
RLP_Mag_Pel_Earned_all{h,1}=RLP_Mag_Pel_Earned;
RLP_Mag_Suc_Earned_all{h,1}=RLP_Mag_Suc_Earned;
RLP_Mag_Pel_Free_all{h,1}=RLP_Mag_Pel_Free;
RLP_Mag_Suc_Free_all{h,1}=RLP_Mag_Suc_Free;

%Extracting signals
LLP_Signal_dFF=[];
LLP_LLP_Signal_dFF=[];
LLP_Mag_Empty_Signal_dFF=[];
LLP_Mag_Pel_Earned_Signal_dFF=[];
LLP_Mag_Suc_Earned_Signal_dFF=[];
LLP_Mag_Pel_Free_Signal_dFF=[];
LLP_Mag_Suc_Free_Signal_dFF=[];

RLP_Signal_dFF=[];
RLP_RLP_Signal_dFF=[];
RLP_Mag_Empty_Signal_dFF=[];
RLP_Mag_Pel_Earned_Signal_dFF=[];
RLP_Mag_Suc_Earned_Signal_dFF=[];
RLP_Mag_Pel_Free_Signal_dFF=[];
RLP_Mag_Suc_Free_Signal_dFF=[];

for k=1:length(Locate_LLP);
   if Locate_LLP(k,1)<= length(Signals_working_dFF(:,1));
   LLP_Signal_dFF(k,:)=Signals_working_dFF(Locate_LLP(k,1),:);    
   end
end

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

for k=1:length(LLP_Mag_Pel_Earned);
    if LLP_Mag_Pel_Earned(k,1)<= length(Signals_working_dFF(:,1));
   LLP_Mag_Pel_Earned_Signal_dFF(k,:)=Signals_working_dFF(LLP_Mag_Pel_Earned(k,1),:);   
    end
end

for k=1:length(LLP_Mag_Suc_Earned);
    if LLP_Mag_Suc_Earned(k,1)<= length(Signals_working_dFF(:,1));
   LLP_Mag_Suc_Earned_Signal_dFF(k,:)=Signals_working_dFF(LLP_Mag_Suc_Earned(k,1),:);   
    end
end

for k=1:length(LLP_Mag_Pel_Free);
    if LLP_Mag_Pel_Free(k,1)<= length(Signals_working_dFF(:,1));
   LLP_Mag_Pel_Free_Signal_dFF(k,:)=Signals_working_dFF(LLP_Mag_Pel_Free(k,1),:);   
    end
end

for k=1:length(LLP_Mag_Suc_Free);
    if LLP_Mag_Suc_Free(k,1)<= length(Signals_working_dFF(:,1));
   LLP_Mag_Suc_Free_Signal_dFF(k,:)=Signals_working_dFF(LLP_Mag_Suc_Free(k,1),:);   
    end
end

for k=1:length(Locate_RLP);
    if Locate_RLP(k,1)<= length(Signals_working_dFF(:,1));
   RLP_Signal_dFF(k,:)=Signals_working_dFF(Locate_RLP(k,1),:);      
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

for k=1:length(RLP_Mag_Pel_Earned);
    if RLP_Mag_Pel_Earned(k,1)<= length(Signals_working_dFF(:,1));
   RLP_Mag_Pel_Earned_Signal_dFF(k,:)=Signals_working_dFF(RLP_Mag_Pel_Earned(k,1),:);     
    end
end

for k=1:length(RLP_Mag_Suc_Earned);
    if RLP_Mag_Suc_Earned(k,1)<= length(Signals_working_dFF(:,1));
   RLP_Mag_Suc_Earned_Signal_dFF(k,:)=Signals_working_dFF(RLP_Mag_Suc_Earned(k,1),:);      
    end
end

for k=1:length(RLP_Mag_Pel_Free);
    if RLP_Mag_Pel_Free(k,1)<= length(Signals_working_dFF(:,1));
   RLP_Mag_Pel_Free_Signal_dFF(k,:)=Signals_working_dFF(RLP_Mag_Pel_Free(k,1),:);     
    end
end

for k=1:length(RLP_Mag_Suc_Free);
    if RLP_Mag_Suc_Free(k,1)<= length(Signals_working_dFF(:,1));
   RLP_Mag_Suc_Free_Signal_dFF(k,:)=Signals_working_dFF(RLP_Mag_Suc_Free(k,1),:);      
    end
end

LLP_Signal_all{h,1}=LLP_Signal_dFF;

LLP_LLP_Signal_all{h,1}=LLP_LLP_Signal_dFF;

LLP_Mag_Empty_Signal_all{h,1}=LLP_Mag_Empty_Signal_dFF;

LLP_Mag_Pel_Earned_Signal_all{h,1}=LLP_Mag_Pel_Earned_Signal_dFF;

LLP_Mag_Suc_Earned_Signal_all{h,1}=LLP_Mag_Suc_Earned_Signal_dFF;

LLP_Mag_Pel_Free_Signal_all{h,1}=LLP_Mag_Pel_Free_Signal_dFF;

LLP_Mag_Suc_Free_Signal_all{h,1}=LLP_Mag_Suc_Free_Signal_dFF;

RLP_Signal_all{h,1}=RLP_Signal_dFF;

RLP_RLP_Signal_all{h,1}=RLP_RLP_Signal_dFF;

RLP_Mag_Empty_Signal_all{h,1}=RLP_Mag_Empty_Signal_dFF;

RLP_Mag_Pel_Earned_Signal_all{h,1}=RLP_Mag_Pel_Earned_Signal_dFF;

RLP_Mag_Suc_Earned_Signal_all{h,1}=RLP_Mag_Suc_Earned_Signal_dFF;

RLP_Mag_Pel_Free_Signal_all{h,1}=RLP_Mag_Pel_Free_Signal_dFF;

RLP_Mag_Suc_Free_Signal_all{h,1}=RLP_Mag_Suc_Free_Signal_dFF;

%///////////////////MAGALIGNED//////////////////////////////
%/////////////////empty mag checks////////////////////////////
magaligned_working=[];
magaligned_LLP_Mag_Empty=[];
for j=1:length(LLP_Mag_Empty);
    magaligned_working(j,1)=find(Locate_Mag_Empty(:,1)>LLP_Mag_Empty(j),1);
    magaligned_LLP_Mag_Empty(j,1)=Locate_Mag_Empty(magaligned_working(j,1),1);
end

magaligned_working=[];
magaligned_RLP_Mag_Empty=[];
for j=1:length(RLP_Mag_Empty);
    magaligned_working(j,1)=find(Locate_Mag_Empty(:,1)>RLP_Mag_Empty(j),1);
    magaligned_RLP_Mag_Empty(j,1)=Locate_Mag_Empty(magaligned_working(j,1),1);
end

magaligned_LLP_Mag_Empty_all{h,1}=magaligned_LLP_Mag_Empty;
magaligned_RLP_Mag_Empty_all{h,1}=magaligned_RLP_Mag_Empty;

magaligned_LLP_Mag_Empty_Signal_dFF=[];
magaligned_RLP_Mag_Empty_Signal_dFF=[];

for k=1:length(magaligned_LLP_Mag_Empty);
    if magaligned_LLP_Mag_Empty(k,1)<= length(Signals_working_dFF(:,1)); 
   magaligned_LLP_Mag_Empty_Signal_dFF(k,:)=Signals_working_dFF(magaligned_LLP_Mag_Empty(k,1),:);    
    end
end

for k=1:length(magaligned_RLP_Mag_Empty);
    if magaligned_RLP_Mag_Empty(k,1)<= length(Signals_working_dFF(:,1)); 
   magaligned_RLP_Mag_Empty_Signal_dFF(k,:)=Signals_working_dFF(magaligned_RLP_Mag_Empty(k,1),:);    
    end
end

magaligned_LLP_Mag_Empty_Signal_all{h,1}=magaligned_LLP_Mag_Empty_Signal_dFF;
magaligned_RLP_Mag_Empty_Signal_all{h,1}=magaligned_RLP_Mag_Empty_Signal_dFF;

%///////////////////Pel/////////////////////////////////
 magaligned_LLP_Mag_Pel_Earned=[];
    for j=1:length(LLP_Mag_Pel_Earned);
    magaligned_LLP_Mag_Pel_Earned(j,1)=First_Mag_Pel_Earned(find(First_Mag_Pel_Earned>LLP_Mag_Pel_Earned(j),1));
    end;

 magaligned_RLP_Mag_Pel_Earned=[];
    for j=1:length(RLP_Mag_Pel_Earned);
    magaligned_RLP_Mag_Pel_Earned(j,1)=First_Mag_Pel_Earned(find(First_Mag_Pel_Earned>RLP_Mag_Pel_Earned(j),1));
    end;

magaligned_LLP_Mag_Pel_Earned_all{h,1}=magaligned_LLP_Mag_Pel_Earned;
magaligned_RLP_Mag_Pel_Earned_all{h,1}=magaligned_RLP_Mag_Pel_Earned;
    
magaligned_LLP_Mag_Pel_Earned_Signal_dFF=[];
magaligned_RLP_Mag_Pel_Earned_Signal_dFF=[];

for k=1:length(magaligned_LLP_Mag_Pel_Earned);
    if magaligned_LLP_Mag_Pel_Earned(k,1)<= length(Signals_working_dFF(:,1)); 
   magaligned_LLP_Mag_Pel_Earned_Signal_dFF(k,:)=Signals_working_dFF(magaligned_LLP_Mag_Pel_Earned(k,1),:);    
    end
end

for k=1:length(magaligned_RLP_Mag_Pel_Earned);
    if magaligned_RLP_Mag_Pel_Earned(k,1)<= length(Signals_working_dFF(:,1)); 
   magaligned_RLP_Mag_Pel_Earned_Signal_dFF(k,:)=Signals_working_dFF(magaligned_RLP_Mag_Pel_Earned(k,1),:);    
    end
end

magaligned_LLP_Mag_Pel_Earned_Signal_all{h,1}=magaligned_LLP_Mag_Pel_Earned_Signal_dFF;
magaligned_RLP_Mag_Pel_Earned_Signal_all{h,1}=magaligned_RLP_Mag_Pel_Earned_Signal_dFF;

%///////////////////Suc/////////////////////////////////
 magaligned_LLP_Mag_Suc_Earned=[];
    for j=1:length(LLP_Mag_Suc_Earned);
    magaligned_LLP_Mag_Suc_Earned(j,1)=First_Mag_Suc_Earned(find(First_Mag_Suc_Earned>LLP_Mag_Suc_Earned(j),1));
    end;

 magaligned_RLP_Mag_Suc_Earned=[];
    for j=1:length(RLP_Mag_Suc_Earned);
    magaligned_RLP_Mag_Suc_Earned(j,1)=First_Mag_Suc_Earned(find(First_Mag_Suc_Earned>RLP_Mag_Suc_Earned(j),1));
    end;

magaligned_LLP_Mag_Suc_Earned_all{h,1}=magaligned_LLP_Mag_Suc_Earned;
magaligned_RLP_Mag_Suc_Earned_all{h,1}=magaligned_RLP_Mag_Suc_Earned;
    
magaligned_LLP_Mag_Suc_Earned_Signal_dFF=[];
magaligned_RLP_Mag_Suc_Earned_Signal_dFF=[];

for k=1:length(magaligned_LLP_Mag_Suc_Earned);
    if magaligned_LLP_Mag_Suc_Earned(k,1)<= length(Signals_working_dFF(:,1)); 
   magaligned_LLP_Mag_Suc_Earned_Signal_dFF(k,:)=Signals_working_dFF(magaligned_LLP_Mag_Suc_Earned(k,1),:);    
    end
end

for k=1:length(magaligned_RLP_Mag_Suc_Earned);
    if magaligned_RLP_Mag_Suc_Earned(k,1)<= length(Signals_working_dFF(:,1)); 
   magaligned_RLP_Mag_Suc_Earned_Signal_dFF(k,:)=Signals_working_dFF(magaligned_RLP_Mag_Suc_Earned(k,1),:);    
    end
end

magaligned_LLP_Mag_Suc_Earned_Signal_all{h,1}=magaligned_LLP_Mag_Suc_Earned_Signal_dFF;
magaligned_RLP_Mag_Suc_Earned_Signal_all{h,1}=magaligned_RLP_Mag_Suc_Earned_Signal_dFF;
    
end

