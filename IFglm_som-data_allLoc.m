load('data\som_lif_glm_data.mat');
load('data\som_lif_glm_voltage_data.mat');

ind_e_all=[3 5];

% ind_sp_all=[1 4]; %cell 1
% ind_sp_all=[18 17]; %cell 7
ind_sp_all=[28 27]; %cell 12

%%

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

locVec_sp=locVec;

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

locVec_sp=locVec;


locVec_sp_pt2=locVec_sp;
locVec_spC_pt2=locVec_spC;

%%
locVec_sp=union(locVec_sp_pt1,locVec_sp_pt2);
locVec_spC=locVec_spC_pt1+locVec_spC_pt2;



%%




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


%%


%% first run R


% betahat_conv_allLoc=[-8.538058 -23.86868 0 0 0 0 0 0 0.004823065 0 0 0 0 0.00313147 0.004746599 0 0 0 0 0.005490498 0 0 0 0 0 0 0];
%%%% cell 1 location 1 all stim (10ms & 4ms) all locations (include zeros)
% betahat_conv_allLoc=[-7.199473 0 0 0 0.001313774 0 0 0 0 0.005691483 0.002042199 0 0 0.001668808 0.00426648 0.007594597 0 0.002106168 0.004907374 0.004581382 0.004905718 0 0 0 0 0 0];    
%%%% cell 7 location 2 all stim all locations (include zeros)
% betahat_conv_allLoc=[-8.578219 -3.394433 0 0 0 0 0 0 0 0.008505008 0 0 0.005525908 0.004630355 0.006062297 0 0 0 0 0 0 0 0 0 0 0 0];
%%%% cell 12 location 1 all stim all locations (include zeros)

gain_vec=zeros(nloc_sp,1);
gain_vec(loc_min)=betahat_conv_allLoc(3);
loc_mmVec=[1:loc_min-1 loc_min+1:nloc_sp];
for i=1:nloc_sp-1
    gain_vec(loc_mmVec(i))=betahat_conv_allLoc(3)+betahat_conv_allLoc(3+i);
end


%% prediction



%% current type 1
ind_sp=ind_sp_all(1);ind_e=ind_e_all(1);

locVec=1:25;Ie_byLoc=[];
for loc=locVec
    Ie=repmat(all_data{ind_e}{loc}(:,1:20:2000)',1,9)*-1;
    Ie_byLoc{loc}=Ie;
end

ntrial=100; perloc=9;
pred_meanSp=zeros(1,nloc);
pred_stdSp=zeros(1,nloc);
pred_latency1st=zeros(1,nloc); %matrix of mean latency for 1st spike: amp*loc
for ii=1:nloc_sp
    loc=locVec_sp(ii);
    Ie_fit=Ie_byLoc{loc}(:,1);
for jj=1
    betahat_conv=zeros(1,3);
    betahat_conv(1:2)=betahat_conv_allLoc(1:2);
    betahat_conv(3)=gain_vec(ii);
    
ntime=length(Ie_fit);
predTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);spiketime=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);step=zeros(ntime,1);
    j=1;lastSpT=0;numSp=0;
    fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
    fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*Ie_fit(lastSpT+1:j);
    
    %%%%link function f(x) = 1+x, x>0; = exp(x), x<0.
    
    step0=betahat_conv(1)+betahat_conv(3)*fit_expg_k(j,tr);
    if step0>0
        step(j)=step0+1;
    else
        step(j)=exp(step0);
    end
    
    lambdaMat(j,tr)=step(j);
    
while (j<=ntime)
    if sum(step)>tao
        step(1:j)=0;
        predTrain(j,tr)=1; %spike
        tao=exprnd(1);
        lastSpT=j;numSp=numSp+1;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*Ie_fit(lastSpT+1:j);
            
            step0=betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr);
            if step0>0
                step(j)=step0+1;
            else
                step(j)=exp(step0);
            end
            
            lambdaMat(j,tr)=step(j);
        end
    else
        predTrain(j,tr)=0;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*Ie_fit(lastSpT+1:j);
            
            if numSp==0
                step0=betahat_conv(1)+betahat_conv(3)*fit_expg_k(j,tr);
                if step0>0
                    step(j)=step0+1;
                else
                    step(j)=exp(step0);
                end
            else
                step0=betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr);
                if step0>0
                    step(j)=step0+1;
                else
                    step(j)=exp(step0);
                end
            end
            
            lambdaMat(j,tr)=step(j);
        end
    end
end
end

[r2,n2]=max(predTrain~=0,[],1);n2(r2==0)=ntime;pred_latency1st(jj,loc)=mean(n2);
pred_meanSp(jj,loc)=round(mean(sum(predTrain,1))*100)/100;
pred_stdSp(jj,loc)=round(std(sum(predTrain,1))*100)/100;

end
predTrain_allLoc{loc}=predTrain;
end
pred_latency1st(pred_latency1st==0)=ntime;

figure;
subplot(241);
    plot(all_data{ind_e}{12}(:,1:20:2000)*-1);
    set(gca,'FontSize',12);
    ylabel({'Current','at cell center'});
    ylim([0 300]);box off;
subplot(242);
    imagesc(reshape(pred_latency1st,5,5)');
    set(gca,'FontSize',12);caxis([0 100]);
    title('Sim: mean 1st spike latency');box off;
    colormap(flipud(hot));%colorbar;
subplot(243);
    imagesc(reshape(pred_meanSp,5,5)');
    set(gca,'FontSize',12);
    caxis([0 4]); %som cell 7,12
%     caxis([0 2]); %som cell 1
    title('Sim: mean');box off;
    colormap(flipud(hot));%colorbar;
subplot(244);
    imagesc(reshape(pred_stdSp,5,5)');
    set(gca,'FontSize',12);caxis([0 2]);
    title('Sim: std dev');box off;
    colormap(flipud(hot));%colorbar;

    
%% current type 2
ind_sp=ind_sp_all(2);ind_e=ind_e_all(2);

locVec=1:25;Ie_byLoc=[];
for loc=locVec
    Ie=repmat(all_data{ind_e}{loc}(:,1:20:2000)',1,9)*-1;
    Ie_byLoc{loc}=Ie;
end

ntrial=100; perloc=9;
pred_meanSp=zeros(1,nloc);
pred_stdSp=zeros(1,nloc);
pred_latency1st=zeros(1,nloc); %matrix of mean latency for 1st spike: amp*loc
for ii=1:nloc_sp
    loc=locVec_sp(ii);
    Ie_fit=Ie_byLoc{loc}(:,1);
for jj=1
    betahat_conv=zeros(1,3);
    betahat_conv(1:2)=betahat_conv_allLoc(1:2);
    betahat_conv(3)=gain_vec(ii);
    
ntime=length(Ie_fit);
predTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);spiketime=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);step=zeros(ntime,1);
    j=1;lastSpT=0;numSp=0;
    fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
    fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*Ie_fit(lastSpT+1:j);
    
    %%%%link function f(x) = 1+x, x>0; = exp(x), x<0.
    
    step0=betahat_conv(1)+betahat_conv(3)*fit_expg_k(j,tr);
    if step0>0
        step(j)=step0+1;
    else
        step(j)=exp(step0);
    end
    
    lambdaMat(j,tr)=step(j);
    
while (j<=ntime)
    if sum(step)>tao
        step(1:j)=0;
        predTrain(j,tr)=1; %spike
        tao=exprnd(1);
        lastSpT=j;numSp=numSp+1;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*Ie_fit(lastSpT+1:j);
            
            step0=betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr);
            if step0>0
                step(j)=step0+1;
            else
                step(j)=exp(step0);
            end
            
            lambdaMat(j,tr)=step(j);
        end
    else
        predTrain(j,tr)=0;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*Ie_fit(lastSpT+1:j);
            
            if numSp==0
                step0=betahat_conv(1)+betahat_conv(3)*fit_expg_k(j,tr);
                if step0>0
                    step(j)=step0+1;
                else
                    step(j)=exp(step0);
                end
            else
                step0=betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr);
                if step0>0
                    step(j)=step0+1;
                else
                    step(j)=exp(step0);
                end
            end
            
            lambdaMat(j,tr)=step(j);
        end
    end
end
end

[r2,n2]=max(predTrain~=0,[],1);n2(r2==0)=ntime;pred_latency1st(jj,loc)=mean(n2);
pred_meanSp(jj,loc)=round(mean(sum(predTrain,1))*100)/100;
pred_stdSp(jj,loc)=round(std(sum(predTrain,1))*100)/100;

end
predTrain_allLoc{loc}=predTrain;
end

pred_latency1st(pred_latency1st==0)=ntime;

subplot(245);
    plot(all_data{ind_e}{12}(:,1:20:2000)*-1);
    set(gca,'FontSize',12);
    ylabel({'Current','at cell center'});
    ylim([0 300]);box off;
subplot(246);
    imagesc(reshape(pred_latency1st,5,5)');
    set(gca,'FontSize',12);caxis([0 100]);
    title('Sim: mean 1st spike latency');box off;
    colormap(flipud(hot));%colorbar;
subplot(247);
    imagesc(reshape(pred_meanSp,5,5)');
    set(gca,'FontSize',12);
    caxis([0 4]); %som cell 7,12 
%     caxis([0 2]); %som cell 1
    title('Sim: mean');box off;
    colormap(flipud(hot));%colorbar;
subplot(248);
    imagesc(reshape(pred_stdSp,5,5)');
    set(gca,'FontSize',12);caxis([0 2]);
    title('Sim: std dev');box off;
    colormap(flipud(hot));%colorbar;

    
%%
figure;
subplot(211);
imagesc(all_data{ind_sp}{12});colormap(flipud(gray));
title('Real data');
xlabel('Time (ms)');
subplot(212);
imagesc(predTrain_allLoc{12}');colormap(flipud(gray));
title('Prediction');
xlabel('Time (ms)');ylabel('Trial');
