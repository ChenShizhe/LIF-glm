load('data\spikes_fix.mat');
load('data\current_template.mat');

%%
%choose a cell
%ce=1; %locVec=[61 50 73 51 63 17 84 66 64];
%ce=2; %locVec=[62 49 39 40 51 61 50 60 29];

ce=6;

locVec=1:121;nloc=length(locVec);

for bl=1:15
    spBlockAll1ms{bl}=spikes_downres{ce}(1+(bl-1)*121:121+(bl-1)*121,:);
end

I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];

%%
real_meanSp=zeros(3,nloc);
real_stdSp=zeros(3,nloc);
real_latency1st=zeros(3,nloc); %matrix of mean latency for 1st spike: amp*loc

spTrain_byLoc=[];locVec_sp=[];
for loc=locVec
    sp_eg_100mV=[spBlockAll1ms{11}(loc,:)' spBlockAll1ms{12}(loc,:)' spBlockAll1ms{13}(loc,:)' spBlockAll1ms{14}(loc,:)' spBlockAll1ms{14}(loc,:)'];
    sp_eg_50mV=[spBlockAll1ms{6}(loc,:)' spBlockAll1ms{7}(loc,:)' spBlockAll1ms{8}(loc,:)' spBlockAll1ms{9}(loc,:)' spBlockAll1ms{10}(loc,:)'];
    sp_eg_25mV=[spBlockAll1ms{1}(loc,:)' spBlockAll1ms{2}(loc,:)' spBlockAll1ms{3}(loc,:)' spBlockAll1ms{4}(loc,:)' spBlockAll1ms{5}(loc,:)'];
    spTrain=[sp_eg_100mV sp_eg_50mV sp_eg_25mV];
    if sum(spTrain(:))>0
        locVec_sp=[locVec_sp;loc];
    end
    spTrain_byLoc{loc}=spTrain;
    
    [r100,n100]=max(sp_eg_100mV~=0,[],1);n100(r100==0)=75;real_latency1st(1,loc)=mean(n100);
    [r50,n50]=max(sp_eg_50mV~=0,[],1);n50(r50==0)=75;real_latency1st(2,loc)=mean(n50);
    [r25,n25]=max(sp_eg_25mV~=0,[],1);n25(r25==0)=75;real_latency1st(3,loc)=mean(n25);
    
    real_meanSp(1,loc)=round(mean(sum(sp_eg_100mV,1))*100)/100;
    real_stdSp(1,loc)=round(std(sum(sp_eg_100mV,1))*100)/100;    
    real_meanSp(2,loc)=round(mean(sum(sp_eg_50mV,1))*100)/100;
    real_stdSp(2,loc)=round(std(sum(sp_eg_50mV,1))*100)/100;    
    real_meanSp(3,loc)=round(mean(sum(sp_eg_25mV,1))*100)/100;
    real_stdSp(3,loc)=round(std(sum(sp_eg_25mV,1))*100)/100;
end

figure;
for pid=1:3
    subplot(1,3,pid);
    imagesc(reshape(real_latency1st(pid,:),11,11)');
    set(gca,'FontSize',16);caxis([0 75]);
    title('Real: mean 1st spike latency');
    colormap(flipud(hot));colorbar;
end

figure;
for pid=1:3
    subplot(1,3,pid);
    imagesc(reshape(real_meanSp(pid,:),11,11)');
    set(gca,'FontSize',16);caxis([0 2]);
    title('Real: spike mean');
    colormap(flipud(hot));colorbar;
end


%%
nloc_sp=length(locVec_sp);

link = @link_test;  %link = @(mu) mu + log(1-exp(-mu));
derlink = @derlink_test;
invlink = @invlink_test;
F = {link, derlink, invlink};

trainM=[];
for loc=locVec_sp
    trainM=[trainM spTrain_byLoc{loc}];
end
I_eg=repmat(I_e,1,nloc_sp);

%% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % with an indicator function, 1 = beyond 1st spike % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%
%%%%with an indicator for >1st spike
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

%%%%double check ind1st is correct:
% figure;
% subplot(121);imagesc(reshape(trainM,75,15*nloc_sp)');colormap(flipud(gray));
% subplot(122);imagesc(ind1st');colormap(flipud(gray));


%% log likelihood
% gVec=0.025:0.025:2;
% logL_vec=zeros(length(gVec),1);
% nii=size(spTrain,1)*size(spTrain,2);
% 
% for i=1:length(gVec)
%     g=gVec(i);
%     [expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time
%     for ii=1:nloc_sp-1
%     indloc=zeros(size(trainM,1)*size(trainM,2),1);
%     indloc((ii-1)*nii+1:ii*nii)=1;
%     expg_k_loc(:,ii)=expg_k(:).*indloc;
%     end
%     [betahat_conv_allLoc,~,stats_conv]=glmfit([ind1st(:).*expg_Vreset(:) expg_k(:) expg_k_loc],trainM(:),'Poisson','link',F);
%     lambdahat_allLoc=glmval(betahat_conv_allLoc,[ind1st(:).*expg_Vreset(:) expg_k(:) expg_k_loc],F);
%     logL=sum(log(poisspdf(trainM(:),lambdahat_allLoc)));
%     logL_vec(i)=logL;
% end
% figure;plot(logL_vec);set(gca,'XTick',5:5:length(gVec),'XTickLabel',gVec(5:5:length(gVec)),'FontSize',16);
% xlabel('g');ylabel('log-likelihood');
% [find(logL_vec==max(logL_vec)) gVec(find(logL_vec==max(logL_vec)))]
% 
%% glm fit
%%%%g_MLE
% g=gVec(find(logL_vec==max(logL_vec)));
g=0.0000001;
[expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time
nii=size(spTrain,1)*size(spTrain,2);
for ii=1:nloc_sp-1
    indloc=zeros(size(trainM,1)*size(trainM,2),1);
    indloc((ii-1)*nii+1:ii*nii)=1;
    expg_k_loc(:,ii)=expg_k(:).*indloc;
end
[betahat_conv_allLoc,~,stats_conv]=glmfit([ind1st(:).*expg_Vreset(:) expg_k(:) expg_k_loc],trainM(:),'Poisson','link',F);
[betahat_conv_allLoc(1); betahat_conv_allLoc(2)]

gain_vec=zeros(length(locVec_sp),1);
gain_vec(end)=betahat_conv_allLoc(3);
for i=1:length(locVec_sp)-1
    gain_vec(i)=betahat_conv_allLoc(3)+betahat_conv_allLoc(3+i);
end
figure;plot(gain_vec,'*');
xlim([0 length(locVec_sp)]);
hold on
plot([0 80],[0 0],'k');
hold off

%% predicted voltage
%%%%for diagnostic purposes, 
%%%%it would be helpful to plot the GLM predicted voltage next to the true spikes,
%%%%especially in cases where the neuron spikes predictably:
%%%%if the predicted voltage consistently gets large in the wrong place,
%%%%this would indicate that the input current model should be changed

% volt_allLoc=[ones(size(trainM(:),1),1) ind1st(:).*expg_Vreset(:) expg_k(:) expg_k_loc]*betahat_conv_allLoc;

% figure;
% subplot(5,1,1:4);
% plot(volt_allLoc(1:7500));
% box off;
% set(gca,'FontSize',16);
% hold on
% plot(trainM(1:7500),'k');
% hold off
% xlim([0 7500]);ylim([-20 2]);
% ylabel('Fitted voltage');
% subplot(5,1,5);
% plot(I_eg(1:7500),'g');box off;
% set(gca,'FontSize',16);
% xlim([0 7500]);
% xlabel('Time (ms)');ylabel('Input current');

% figure;%suptitle('Estimated voltage (by location)');
% for ii=1:nloc_sp
%     subplottight(11,11,locVec_sp(ii));
%     axis off;
%     volt_byLoc=volt_allLoc(75*15*(ii-1)+1:75*15*ii);
%     sp_byLoc=trainM(75*15*(ii-1)+1:75*15*ii);
%     volt=reshape(volt_byLoc,[75 15]);
%     sp=reshape(sp_byLoc,[75 15]);
%     hold on;
%     imagesc(1:75,14.5:-1:0.5,sp');colormap(flipud(gray));
%     set(gca,'ydir','normal');
%     for jj=1:15
%         plot(volt(:,jj)-min(volt(:,jj))+(14-jj),'b');
%     end
%     hold off
%     xlim([0 75]);ylim([-1 18]);
% end

%% prediction statistics: first spike latency; mean number of spike

ntrial=100;
pred_meanSp=zeros(3,nloc);
pred_stdSp=zeros(3,nloc);
pred_latency1st=zeros(3,nloc); %matrix of mean latency for 1st spike: amp*loc
for ii=1:nloc_sp
    loc=locVec_sp(ii);
for jj=1:3
    betahat_conv=zeros(1,3);
    betahat_conv(1:2)=betahat_conv_allLoc(1:2);
    betahat_conv(3)=gain_vec(ii);
    
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

[r2,n2]=max(predTrain~=0,[],1);n2(r2==0)=75;pred_latency1st(jj,loc)=mean(n2);
pred_meanSp(jj,loc)=round(mean(sum(predTrain,1))*100)/100;
pred_stdSp(jj,loc)=round(std(sum(predTrain,1))*100)/100;
end
end
pred_latency1st(pred_latency1st==0)=75;

figure;
for pid=1:3
    subplot(1,3,pid);
    imagesc(reshape(pred_latency1st(pid,:),11,11)');
    set(gca,'FontSize',16);caxis([0 75]);
    title('Sim: mean 1st spike latency');
    colormap(flipud(hot));colorbar;
end

figure;
for pid=1:3
    subplot(1,3,pid);
    imagesc(reshape(pred_meanSp(pid,:),11,11)');
    set(gca,'FontSize',16);caxis([0 2]);
    title('Sim: spike mean');
    colormap(flipud(hot));colorbar;
end
