load('data\spikes_fix.mat');
load('data\current_template.mat');

%%
%choose a cell
%ce=1; locVec=[61 50 73 51 63 17 84 66 64];
ce=2; locVec=[62 49 39 40 51 61 50 60 29];

nloc=length(locVec);

for bl=1:15
    spBlockAll1ms{bl}=spikes_downres{ce}(1+(bl-1)*121:121+(bl-1)*121,:);
end

I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];

spTrain_allLoc=[];
Ie_allLoc=[];
ind_allLoc=zeros(size(I_e,1)*size(I_e,2),9);


real_meanSp=zeros(3,9);
real_stdSp=zeros(3,9);
real_latency1st=zeros(3,9); %matrix of mean latency for 1st spike: amp*loc

spTrain_byLoc=[];
figure;suptitle('Real data (by location)');
for ii=1:nloc
    loc=locVec(ii);
    sp_eg_100mV=[spBlockAll1ms{11}(loc,:)' spBlockAll1ms{12}(loc,:)' spBlockAll1ms{13}(loc,:)' spBlockAll1ms{14}(loc,:)' spBlockAll1ms{14}(loc,:)'];
    sp_eg_50mV=[spBlockAll1ms{6}(loc,:)' spBlockAll1ms{7}(loc,:)' spBlockAll1ms{8}(loc,:)' spBlockAll1ms{9}(loc,:)' spBlockAll1ms{10}(loc,:)'];
    sp_eg_25mV=[spBlockAll1ms{1}(loc,:)' spBlockAll1ms{2}(loc,:)' spBlockAll1ms{3}(loc,:)' spBlockAll1ms{4}(loc,:)' spBlockAll1ms{5}(loc,:)'];
    spTrain=[sp_eg_100mV sp_eg_50mV sp_eg_25mV];
    spTrain_byLoc{loc}=spTrain;
    
    [r100,n100]=max(sp_eg_100mV~=0,[],1);real_latency1st(1,ii)=mean(n100);
    [r50,n50]=max(sp_eg_50mV~=0,[],1);real_latency1st(2,ii)=mean(n50);
    [r25,n25]=max(sp_eg_25mV~=0,[],1);real_latency1st(3,ii)=mean(n25);
    
    real_meanSp(1,ii)=round(mean(sum(sp_eg_100mV,1))*100)/100;
    real_stdSp(1,ii)=round(std(sum(sp_eg_100mV,1))*100)/100;    
    real_meanSp(2,ii)=round(mean(sum(sp_eg_50mV,1))*100)/100;
    real_stdSp(2,ii)=round(std(sum(sp_eg_50mV,1))*100)/100;    
    real_meanSp(3,ii)=round(mean(sum(sp_eg_25mV,1))*100)/100;
    real_stdSp(3,ii)=round(std(sum(sp_eg_25mV,1))*100)/100;    

    spTrain_allLoc=[spTrain_allLoc;spTrain(:)];
    Ie_allLoc=[Ie_allLoc;I_e(:)];
    nii=size(spTrain,1)*size(spTrain,2);
    ind_allLoc((ii-1)*nii+1:ii*nii,ii)=1;
    
    subplot(3,3,ii);
    imagesc(spTrain');
    title(['Loc=' num2str(locVec(ii))]);
    hold on
    plot([5 5],[0 122],'k:');
    plot([15 15],[0 122],'k:');
    plot([0 75],[5.5 5.5],'k-');
    plot([0 75],[10.5 10.5],'k-');
    plot([0 75],[15.5 15.5],'k-');
    hold off
    if (ii==1 | ii==4 | ii==7)
        ylabel('25mV 50mV 100mV');
    end
    set(gca,'FontSize',12);
    colormap(flipud(gray));
    xlabel('Time (ms)');
end

real_latency1st
real_meanSp
real_stdSp


spTrain_byStim=[];
for jj=1:15
    for ii=1:9
        loc=locVec(ii);
        spTrain_byStim{jj}(:,ii)=spTrain_byLoc{loc}(:,jj);
    end
end
figure;suptitle('Real data (by amplitude)');
for jj=1:15
    subplot(3,5,jj);
    imagesc(spTrain_byStim{jj}');
    set(gca,'FontSize',12,'YTick',1:1:9,'YTickLabel',locVec);
    hold on
    plot([5 5],[0 122],'k:');
    plot([15 15],[0 122],'k:');
    hold off
    if jj==1
        ylabel('Stim=100mV');
    elseif jj==6
        ylabel('Stim=50mV');
    elseif jj==11
        ylabel('Stim=25mV');
    end
    colormap(flipud(gray));
    xlabel('Time (ms)');
end


%%
link = @link_test;  %link = @(mu) mu + log(1-exp(-mu));
derlink = @derlink_test;
invlink = @invlink_test;
F = {link, derlink, invlink};

trainM=[];
for ii=1:nloc
    loc=locVec(ii);
    trainM=[trainM spTrain_byLoc{loc}];
end
I_eg=repmat(I_e,1,nloc);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % with an indicator function, 1 = beyond 1st spike % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%
%%%%with an indicator for >1st spike
ind1st_allLoc=[];
for ii=1:nloc
    loc=locVec(ii);
    spTrain=spTrain_byLoc{loc};
    ind1st=zeros(size(spTrain));
    for j=1:size(spTrain,2)
        sp1st=find(spTrain(:,j),1);
        ind1st(sp1st+1:end,j)=1; %%%%indicator for >1st spike
    end
    ind1st_allLoc{ii}=ind1st;
end
ind1st=[];
for ii=1:nloc
    ind1st=[ind1st ind1st_allLoc{ii}];
end

%%%%double check ind1st is correct:
figure;
subplot(121);imagesc(reshape(spTrain_allLoc,75,135)');colormap(flipud(gray));
subplot(122);imagesc(ind1st');colormap(flipud(gray));

%% log likelihood
gVec=0.02:0.02:2;
logL_vec=zeros(length(gVec),1);
for i=1:length(gVec)
    g=gVec(i);
    [expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time
    for iloc=1:nloc-1
        expg_k_loc(:,iloc)=expg_k(:).*ind_allLoc(:,iloc+1);
    end
    [betahat_conv_allLoc,~,stats_conv]=glmfit([ind1st(:).*expg_Vreset(:) expg_k(:) expg_k_loc],spTrain_allLoc,'Poisson','link',F);
    lambdahat_allLoc=glmval(betahat_conv_allLoc,[ind1st(:).*expg_Vreset(:) expg_k(:) expg_k_loc],F);
    logL=sum(log(poisspdf(trainM(:),lambdahat_allLoc)));
    logL_vec(i)=logL;
end
figure;plot(logL_vec);set(gca,'XTick',1:5:length(gVec),'XTickLabel',gVec(1:5:end),'FontSize',16);
xlabel('g');ylabel('log-likelihood');
[find(logL_vec==max(logL_vec)) gVec(find(logL_vec==max(logL_vec)))]

%% glm fit
%%%%g_MLE
g=gVec(find(logL_vec==max(logL_vec)));
%g=0.1;
[expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time
for iloc=1:nloc-1
    expg_k_loc(:,iloc)=expg_k(:).*ind_allLoc(:,iloc+1);
end
[betahat_conv_allLoc,~,stats_conv]=glmfit([ind1st(:).*expg_Vreset(:) expg_k(:) expg_k_loc],spTrain_allLoc,'Poisson','link',F);
[betahat_conv_allLoc(1); betahat_conv_allLoc(2);betahat_conv_allLoc(3);betahat_conv_allLoc(3)+betahat_conv_allLoc(4);betahat_conv_allLoc(3)+betahat_conv_allLoc(5);betahat_conv_allLoc(3)+betahat_conv_allLoc(6);betahat_conv_allLoc(3)+betahat_conv_allLoc(7);betahat_conv_allLoc(3)+betahat_conv_allLoc(8);betahat_conv_allLoc(3)+betahat_conv_allLoc(9);betahat_conv_allLoc(3)+betahat_conv_allLoc(10);betahat_conv_allLoc(3)+betahat_conv_allLoc(11)]'

%% predicted voltage
%%%%for diagnostic purposes, 
%%%%it would be helpful to plot the GLM predicted voltage next to the true spikes,
%%%%especially in cases where the neuron spikes predictably:
%%%%if the predicted voltage consistently gets large in the wrong place,
%%%%this would indicate that the input current model should be changed

lambdahat_allLoc=glmval(betahat_conv_allLoc,[ind1st(:).*expg_Vreset(:) expg_k(:) expg_k_loc],F);
figure;suptitle('Estimated voltage (by location)');
for ii=1:nloc
    subplot(3,3,ii);
    lambdahat_byLoc=lambdahat_allLoc(75*15*(ii-1)+1:75*15*ii);
    lambdahat=reshape(lambdahat_byLoc,[75 15]);
    hold on;
    for jj=1:15
        plot(lambdahat(:,jj)+(15-jj)*0.05,'b');
    end
    title(['Loc=' num2str(locVec(ii))]);
    hold on
    plot([5 5],[0 122],'k:');
    plot([15 15],[0 122],'k:');
    plot([0 75],[0.25 0.25],'k');
    plot([0 75],[0.5 0.5],'k');
    hold off
    xlim([0 75]);ylim([0 0.75]);
    if (ii==1 | ii==4 | ii==7)
        ylabel('25mV 50mV 100mV');
    end
end





%% prediction
figure;suptitle('Prediction (by amplitude)');
ntrial=1;
for jj=1:15
    simTrainMat=zeros(75,9);
for ii=1:9
    betahat_conv=zeros(1,3);
    betahat_conv(1:2)=betahat_conv_allLoc(1:2);
    if ii==1
        betahat_conv(3)=betahat_conv_allLoc(3);
    else
        betahat_conv(3)=betahat_conv_allLoc(3)+betahat_conv_allLoc(2+ii);
    end

%% prediction
I_eg_fit=I_e(:,jj);
ntime=length(I_eg_fit);
predTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);spiketime=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);step=zeros(ntime,1);
    j=1;lastSpT=0;numSp=0;
    fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
    fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
    
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
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
            
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
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
            
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

simTrainMat(:,ii)=predTrain;
end
simTrainAll{jj}=simTrainMat;
subplot(3,5,jj);
imagesc(simTrainMat');
if jj==1
    ylabel('Stim=100mV');
elseif jj==6
    ylabel('Stim=50mV');
elseif jj==11
    ylabel('Stim=25mV');
end
set(gca,'YTick',1:1:9,'YTickLabel',locVec,'XTick',0:25:75,'XTickLabel',0:25:75,'FontSize',12);
hold on
plot([5 5],[0 122],'k:');
plot([15 15],[0 122],'k:');
hold off
colormap(flipud(gray));
xlabel('Time (ms)');
box off;
end

simTrain_byLoc=[];
for jj=1:15
    for ii=1:9
    simTrain_byLoc{ii}(:,jj)=simTrainAll{jj}(:,ii);
    end
end
figure;suptitle('Prediction (by location)');
for ii=1:9
    subplot(3,3,ii);
    imagesc(simTrain_byLoc{ii}');
    title(['Loc=' num2str(locVec(ii))]);
    hold on
    plot([5 5],[0 122],'k:');
    plot([15 15],[0 122],'k:');
    plot([0 75],[5.5 5.5],'k-');
    plot([0 75],[10.5 10.5],'k-');
    plot([0 75],[15.5 15.5],'k-');
    hold off
    if (ii==1 | ii==4 | ii==7)
        ylabel('25mV 50mV 100mV');
    end
    set(gca,'FontSize',12);
    colormap(flipud(gray));
    xlabel('Time (ms)');
end


%% prediction statistics: first spike latency; mean number of spike

ntrial=1000;
pred_meanSp=zeros(3,9);
sim_stdSp=zeros(3,9);
pred_latency1st=zeros(3,9); %matrix of mean latency for 1st spike: amp*loc
for ii=1:9
for jj=1:3
    betahat_conv=zeros(1,3);
    betahat_conv(1:2)=betahat_conv_allLoc(1:2);
    if ii==1
        betahat_conv(3)=betahat_conv_allLoc(3);
    else
        betahat_conv(3)=betahat_conv_allLoc(3)+betahat_conv_allLoc(2+ii);
    end

%% prediction
I_eg_fit=I_e(:,(jj-1)*5+1);
ntime=length(I_eg_fit);
predTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);spiketime=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);step=zeros(ntime,1);
    j=1;lastSpT=0;numSp=0;
    fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
    fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
    
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
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
            
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
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
            
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

[r2,n2]=max(predTrain~=0,[],1);
pred_latency1st(jj,ii)=mean(n2);
pred_meanSp(jj,ii)=round(mean(sum(predTrain,1))*100)/100;
pred_stdSp(jj,ii)=round(std(sum(predTrain,1))*100)/100;
end
end

pred_latency1st
pred_meanSp
pred_stdSp

