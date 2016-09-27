load('data\spikes_fix.mat');
load('data\current_template.mat');

%%
%choose a cell
%ce=1; locVec=[61 50 73 51 63 17 84 66 64];
ce=2; locVec=[62 49 39 40 51 61 50 60 29];

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
for ii=1:9
    loc=locVec(ii);
    sp_eg_100mV=[spBlockAll1ms{11}(loc,:)' spBlockAll1ms{12}(loc,:)' spBlockAll1ms{13}(loc,:)' spBlockAll1ms{14}(loc,:)' spBlockAll1ms{14}(loc,:)'];
    sp_eg_50mV=[spBlockAll1ms{6}(loc,:)' spBlockAll1ms{7}(loc,:)' spBlockAll1ms{8}(loc,:)' spBlockAll1ms{9}(loc,:)' spBlockAll1ms{10}(loc,:)'];
    sp_eg_25mV=[spBlockAll1ms{1}(loc,:)' spBlockAll1ms{2}(loc,:)' spBlockAll1ms{3}(loc,:)' spBlockAll1ms{4}(loc,:)' spBlockAll1ms{5}(loc,:)'];
    spTrain=[sp_eg_100mV sp_eg_50mV sp_eg_25mV];
    spTrain_byLoc{loc}=spTrain;
    
    spTrain_allLoc=[spTrain_allLoc;spTrain(:)];
    Ie_allLoc=[Ie_allLoc;I_e(:)];
    nii=size(spTrain,1)*size(spTrain,2);
    ind_allLoc((ii-1)*nii+1:ii*nii,ii)=ii-1;
end

%%
link = @(mu) log(exp(mu)-1);
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

trainM=[];
for ii=1:9
    loc=locVec(ii);
    trainM=[trainM spTrain_byLoc{loc}];
end
I_eg=repmat(I_e,1,9);
g=0.6;
[expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time

%%%%all loc with an indicator function
[betahat_conv_allLoc,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:) ind_allLoc(:,2:9)],spTrain_allLoc,'Poisson','link',F);
betahat_conv_allLoc'


%% prediction
link = @(mu) log(exp(mu)-1);
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

figure;suptitle('Prediction');
simTrainMat=zeros(75,9);
for jj=1:15
for ii=1:9
    betahat_conv=zeros(1,3);
    betahat_conv(1:3)=betahat_conv_allLoc(1:3);
    if jj==1
        betahat_conv(4)=0;
    else
        betahat_conv(4)=betahat_conv_allLoc(ii+2);
    end

%% prediction
I_eg_fit=I_e(:,jj);
ntrial=1;
ntime=length(I_eg_fit);
simTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);spiketime=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);
    j=1;lastSpT=0;
    fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
    fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
    step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr)+betahat_conv(4))+1);
    %step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_k(j,tr)+betahat_conv(3))+1);
    lambdaMat(j,tr)=step(j);
while (j<=ntime)
    if sum(step)>tao
        step(1:j)=0;
        simTrain(j,tr)=1; %spike
        tao=exprnd(1);
        lastSpT=j;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
            step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr)+betahat_conv(4))+1);
            %step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_k(j,tr)+betahat_conv(3))+1);
            lambdaMat(j,tr)=step(j);
        end
    else
        simTrain(j,tr)=0;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
            step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr)+betahat_conv(4))+1);
            %step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_k(j,tr)+betahat_conv(3))+1);
            lambdaMat(j,tr)=step(j);
        end
    end
end
clear tao step
end

simTrainMat(:,ii)=simTrain;
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
set(gca,'YTick',1:1:9,'YTickLabel',locVec,'XTick',0:25:75,'XTickLabel',0:25:75);
hold on
plot([5 5],[0 122],'k:');
plot([15 15],[0 122],'k:');
hold off
set(gca,'FontSize',12);
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
figure;suptitle('Prediction');
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

%%

gVec=0.02:0.02:2;
logL_vec=zeros(length(gVec),1);
for i=1:length(gVec)
    g=gVec(i);
    [expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time
    [betahat_conv_allLoc,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:) ind_allLoc(:,2:9)],spTrain_allLoc,'Poisson','link',F);
    lambdahat_allLoc=glmval(betahat_conv_allLoc,[expg_Vreset(:) expg_k(:) ind_allLoc(:,2:9)],F);
    logL=sum(log(poisspdf(trainM(:),lambdahat_allLoc)));
    logL_vec(i)=logL;
end
figure;plot(logL_vec);set(gca,'XTick',1:5:length(gVec),'XTickLabel',gVec(1:5:end),'FontSize',16);
xlabel('g');ylabel('log-likelihood');
[find(logL_vec==max(logL_vec)) gVec(find(logL_vec==max(logL_vec)))]

