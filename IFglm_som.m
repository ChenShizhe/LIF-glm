function [ design1,design2,Y,nloc_sp,loc_min ] = IFglm_som( ind_sp_all,ind_e_all, loc_sp, e_max, all_data )
%LIF-glmnet for som data, matlab portion
%%% input:
%%%1. ind_sp_all: which spike train, per Excel id;
%%%2. ind_e_all: which current input, per Excel id.
%%%3. loc_sp: if 1, use spiking loc only; if 0, use all loc including zeros. 
%%%4. e_max: if 1, use input current at max loc only; if 0, use all.
%%%5. all_data: som data set.


%% select a cell and an input current
ind_sp=ind_sp_all(1);ind_e=ind_e_all(1);

locVec=1:25;nloc=length(locVec);

locVec_sp=[];locVec_spC=zeros(nloc,1);
for loc=locVec
    spTrain=all_data{ind_sp}{loc}';
    if sum(spTrain(:))>0
        locVec_sp=[locVec_sp;loc];
        locVec_spC(loc)=sum(spTrain(:));
    end
end

if loc_sp==0
    locVec_sp=locVec;
elseif loc_sp==1
    locVec_sp=locVec_sp;
end

locVec_sp_pt1=locVec_sp;
locVec_spC_pt1=locVec_spC;


%% select a cell and an input current
ind_sp=ind_sp_all(2);ind_e=ind_e_all(2);

locVec=1:25;nloc=length(locVec);

locVec_sp=[];locVec_spC=zeros(nloc,1);
for loc=locVec
    spTrain=all_data{ind_sp}{loc}';    
    if sum(spTrain(:))>0
        locVec_sp=[locVec_sp;loc];
        locVec_spC(loc)=sum(spTrain(:));
    end
end

if loc_sp==0
    locVec_sp=locVec;
elseif loc_sp==1
    locVec_sp=locVec_sp;
end

locVec_sp_pt2=locVec_sp;
locVec_spC_pt2=locVec_spC;

%%
locVec_sp=union(locVec_sp_pt1,locVec_sp_pt2);
locVec_spC=locVec_spC_pt1+locVec_spC_pt2;


%% glmnet prep

%% first current type
ind_sp=ind_sp_all(1);ind_e=ind_e_all(1);

locVec=1:25;nloc=length(locVec);

spTrain=[];spTrain_byLoc=[];
Ie_byLoc=[];
for loc=locVec
    spTrain=all_data{ind_sp}{loc}';
    Ie=repmat(all_data{ind_e}{loc}(:,1:20:2000)',1,9)*-1;
    
    spTrain_byLoc{loc}=spTrain;
    Ie_byLoc{loc}=Ie;
    
    [rVol,nVol]=max(spTrain~=0,[],1);nVol(rVol==0)=100;real_latency1st(loc)=mean(nVol);
    real_meanSp(loc)=round(mean(sum(spTrain,1))*100)/100;
    real_stdSp(loc)=round(std(sum(spTrain,1))*100)/100;    
end


nloc_sp=length(locVec_sp);

trainM=[];IeM=[];
for loc=locVec_sp
    trainM=[trainM spTrain_byLoc{loc}];
    IeM=[IeM Ie_byLoc{loc}]; % current known at every location
end

if e_max==1
    IeM=repmat(Ie_byLoc{12},1,nloc_sp); % fixed current at max
elseif e_max==0
    IeM=IeM;
end

% with an indicator for >1st spike
ind1st_allLoc=[];ind1st=[];
for iloc=1:nloc_sp
    loc=locVec_sp(iloc);
    spTrain=spTrain_byLoc{loc};
    ind1st0=zeros(size(spTrain));
    for j=1:size(spTrain,2)
        sp1st=find(spTrain(:,j),1);
        ind1st0(sp1st+1:end,j)=1; %%%%indicator for >1st spike
    end
    ind1st_allLoc{loc}=ind1st0;
    ind1st=[ind1st ind1st0];
end

% glm fit prep
g=0.1;
nii=size(spTrain,1)*size(spTrain,2);
[expg_Vreset,expg_EL,expg_k]=gconv(IeM,trainM,g); %temporally convolve paramters with g upto spike time

% for R
expg_k_loc0=zeros(size(trainM(:),1),nloc_sp);
for ii=1:nloc_sp
    indloc=zeros(size(trainM,1)*size(trainM,2),1);
    indloc((ii-1)*nii+1:ii*nii)=1;
    expg_k_loc0(:,ii)=expg_k(:).*indloc;
end
% loc_min=find(locVec_sp==find(locVec_spC==min(locVec_spC(find(locVec_spC))),1));
loc_min=1;
expg_k_loc=expg_k_loc0(:,[1:loc_min-1 loc_min+1:end]);

design1=[ones(size(expg_k_loc,1),1) ind1st(:).*expg_Vreset(:)];
design2=[expg_k(:) expg_k_loc];
Y=trainM(:);

design1_pt1=design1;
design2_pt1=design2;
Y_pt1=Y;

%%


%% second current type
ind_sp=ind_sp_all(2);ind_e=ind_e_all(2);

locVec=1:25;nloc=length(locVec);

spTrain=[];spTrain_byLoc=[];
Ie_byLoc=[];
for loc=locVec
    spTrain=all_data{ind_sp}{loc}';
    Ie=repmat(all_data{ind_e}{loc}(:,1:20:2000)',1,9)*-1;
    
    spTrain_byLoc{loc}=spTrain;
    Ie_byLoc{loc}=Ie;
    
    [rVol,nVol]=max(spTrain~=0,[],1);nVol(rVol==0)=100;real_latency1st(loc)=mean(nVol);
    real_meanSp(loc)=round(mean(sum(spTrain,1))*100)/100;
    real_stdSp(loc)=round(std(sum(spTrain,1))*100)/100;    
end

trainM=[];IeM=[];
for loc=locVec_sp
    trainM=[trainM spTrain_byLoc{loc}];
    IeM=[IeM Ie_byLoc{loc}]; % current known at every location
end

if e_max==1
    IeM=repmat(Ie_byLoc{12},1,nloc_sp); % fixed current at max
elseif e_max==0
    IeM=IeM;
end

% with an indicator for >1st spike
ind1st_allLoc=[];ind1st=[];
for iloc=1:nloc_sp
    loc=locVec_sp(iloc);
    spTrain=spTrain_byLoc{loc};
    ind1st0=zeros(size(spTrain));
    for j=1:size(spTrain,2)
        sp1st=find(spTrain(:,j),1);
        ind1st0(sp1st+1:end,j)=1; %%%%indicator for >1st spike
    end
    ind1st_allLoc{loc}=ind1st0;
    ind1st=[ind1st ind1st0];
end

% glm fit prep
g=0.1;
nii=size(spTrain,1)*size(spTrain,2);
[expg_Vreset,expg_EL,expg_k]=gconv(IeM,trainM,g); %temporally convolve paramters with g upto spike time

% for R
expg_k_loc0=zeros(size(trainM(:),1),nloc_sp);
for ii=1:nloc_sp
    indloc=zeros(size(trainM,1)*size(trainM,2),1);
    indloc((ii-1)*nii+1:ii*nii)=1;
    expg_k_loc0(:,ii)=expg_k(:).*indloc;
end
% loc_min=find(locVec_sp==find(locVec_spC==min(locVec_spC(find(locVec_spC))),1));
loc_min=1;
expg_k_loc=expg_k_loc0(:,[1:loc_min-1 loc_min+1:end]);

design1=[ones(size(expg_k_loc,1),1) ind1st(:).*expg_Vreset(:)];
design2=[expg_k(:) expg_k_loc];
Y=trainM(:);

design1_pt2=design1;
design2_pt2=design2;
Y_pt2=Y;

%%
clear design1 design2 Y
design1=[design1_pt1; design1_pt2];
design2=[design2_pt1; design2_pt2];
Y=[Y_pt1;Y_pt2];

end

