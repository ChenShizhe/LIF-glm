load('data\current_template.mat');
I_eg_100mV=abs(norm_average_current(1:20:20*75).*100);
I_eg_50mV=abs(norm_average_current(1:20:20*75).*50);
I_eg_25mV=abs(norm_average_current(1:20:20*75).*25);
I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];

% figure;
% plot(I_e);
% set(gca,'FontSize',12);
% xlim([0 75]);
% xlabel('Time (ms)');ylabel('Input current (mV)');
% legend('100 mV','50 mV','25 mV','location','northeast');

%% simulation

figure;suptitle('Simulation (by location)');
for locid=1:9
dt=1; %time step ms
t_end=75;t_StimStart=5;t_StimEnd=15;
V_th=-25; %spike threshold [mV]
E_L=-28; %resting membrane potential [mV]
V_reset=-80; %value to reset voltage to after a spike [mV]
g=0.1; %membrane time constant [ms]
k=0.04-0.002*locid; %membrane resistance [MOhm]
k_vec(locid)=k;

stoc_mu=0;stoc_sigma=0.3;

%%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect=0:dt:t_end; 
V_vect=zeros(1,length(t_vect));V_plot_vect=zeros(1,length(t_vect));
PlotNum=0;
I_Stim_vect=[100 100 100 100 100 50 50 50 50 50 25 25 25 25 25];
spTrain=zeros(t_end,length(I_Stim_vect));

for I_Stim=I_Stim_vect; %loop over different I_Stim values
    PlotNum=PlotNum+1;
    i=1; %index denoting which element of V is being assigned
    
    V_vect(i)=E_L; %first element of V, i.e. value of V at t=0
    
    %%%%square current
%     I_e_vect=zeros(1,t_StimStart/dt); %portion of I_e_vect from t=0 to t_StimStart
%     I_e_vect=[I_e_vect I_Stim*ones(1,1+((t_StimEnd-t_StimStart)/dt))];
%     I_e_vect=[I_e_vect zeros(1,(t_end-t_StimEnd)/dt)];

    %%%%chi-sq shape current
    I_e_vect=[0;I_e(:,PlotNum)];
    
    I_e_vect_mat(:,PlotNum)=I_e_vect;
    
    NumSpikes=0; %holds number of spikes that have occurred
    for t=dt:dt:t_end %loop through values of t in steps of df ms        
%         V_inf = E_L + I_e_vect(i)*R_m; %value that V_vect is exponentially
%                                          %decaying towards at this time step
%         V_vect(i+1) = V_inf + (V_vect(i)-V_inf)*exp(-dt/tau);
%         V_vect(i+1) = V_vect(i) + (E_L-V_vect(i) + I_e_vect(i)*R_m)/tau*dt + sqrt(dt)*normrnd(stoc_mu,stoc_sigma);
        
        V_vect(i+1) = V_vect(i) + ((E_L-V_vect(i))*g + I_e_vect(i)*k)*dt + sqrt(dt)*normrnd(stoc_mu,stoc_sigma); 
        
        %if statement below says what to do if voltage crosses threshold
        if (V_vect(i+1)>V_th) %cell spiked
            V_vect(i+1)=V_reset; %set voltage back to V_reset
            NumSpikes=NumSpikes+1; %add 1 to the total spike count
            spTrain(i,PlotNum)=1;
        end
        i=i+1; %add 1 to index,corresponding to moving forward 1 time step
    end    
end

spTrain_byLoc{locid}=spTrain;
Ie_byLoc{locid}=I_e_vect_mat(1:end-1,:);

subplot(3,3,locid)
imagesc(spTrain');colormap(flipud(gray));
set(gca,'YTick',5:5:15,'YTickLabel',5:5:15,'XTick',0:25:75,'XTickLabel',0:25:75,'FontSize',12);
hold on
plot([t_StimStart t_StimStart],[0 20],'k:');
plot([t_StimEnd t_StimEnd],[0 20],'k:');
plot([0 75],[5.5 5.5],'k-');
plot([0 75],[10.5 10.5],'k-');
plot([0 75],[15.5 15.5],'k-');
hold off
xlabel('Time (ms)');
if (locid==1 | locid==4 | locid==7)
    ylabel('25mV 50mV 100mV');
end

[r100,n100]=max(spTrain(:,1:5)~=0,[],1);real_latency1st(1,locid)=mean(n100.*r100);
[r50,n50]=max(spTrain(:,6:10)~=0,[],1);real_latency1st(2,locid)=mean(n50.*r50);
[r25,n25]=max(spTrain(:,11:15)~=0,[],1);real_latency1st(3,locid)=mean(n25.*r25);
    
real_meanSp(1,locid)=round(mean(sum(spTrain(:,1:5),1))*100)/100;
real_stdSp(1,locid)=round(std(sum(spTrain(:,1:5),1))*100)/100;  
real_meanSp(2,locid)=round(mean(sum(spTrain(:,5:10),1))*100)/100;
real_stdSp(2,locid)=round(std(sum(spTrain(:,5:10),1))*100)/100;    
real_meanSp(3,locid)=round(mean(sum(spTrain(:,11:15),1))*100)/100;
real_stdSp(3,locid)=round(std(sum(spTrain(:,11:15),1))*100)/100;
end

simTrain_byStim=[];
for jj=1:15
    for ii=1:9
    simTrain_byStim{jj}(:,ii)=spTrain_byLoc{ii}(:,jj);
    end
end
figure;suptitle('Stimulation (by amplitude)');
for ii=1:15
    subplot(3,5,ii);
    imagesc(simTrain_byStim{ii}');
    set(gca,'FontSize',12,'YTick',3:3:9,'YTickLabel',3:3:9);
    hold on
    plot([5 5],[0 122],'k:');
    plot([15 15],[0 122],'k:');
    hold off
    if ii==1
        ylabel('Stim=100mV');
    elseif ii==6
        ylabel('Stim=50mV');
    elseif ii==11
        ylabel('Stim=25mV');
    end
    colormap(flipud(gray));
    xlabel('Time (ms)');
end

%%
figure;
for pid=1:3
    subplot(1,3,pid);
    imagesc(reshape(real_latency1st(pid,:),3,3)');
    set(gca,'FontSize',16);caxis([10 20]);
    title('Real: mean 1st spike latency');
    colormap(jet);colorbar;
end

figure;
for pid=1:3
    subplot(1,3,pid);
    imagesc(reshape(real_meanSp(pid,:),3,3)');
    set(gca,'FontSize',16);caxis([0 4]);
    title('Real: spike mean');
    colormap(flipud(hot));colorbar;
end
figure;
for pid=1:3
    subplot(1,3,pid);
    imagesc(reshape(real_stdSp(pid,:),3,3)');
    set(gca,'FontSize',16);caxis([0 1]);
    title('Real: spike std. dev.');
    colormap(flipud(hot));colorbar;
end


%% glm fit prep

%%%%%location dummy variable ind_allLoc
spTrain_allLoc=[];
Ie_allLoc=[];
ind1st_allLoc=[];
ind_allLoc=zeros(size(spTrain,1)*size(spTrain,2),9);
for ii=1:9
    spTrain=spTrain_byLoc{ii};
    ind1st=zeros(size(spTrain));
    for j=1:size(spTrain,2)
        sp1st=find(spTrain(:,j),1);
        ind1st(sp1st+1:end,j)=1; %%%%indicator for >1st spike
    end
    ind1st_allLoc{ii}=ind1st;
    spTrain_allLoc=[spTrain_allLoc;spTrain(:)];
    I_e=Ie_byLoc{ii};
    Ie_allLoc=[Ie_allLoc;I_e(:)];
    nii=size(spTrain,1)*size(spTrain,2);
    ind_allLoc((ii-1)*nii+1:ii*nii,ii)=1;
end

% link = @(mu) log(exp(mu)-1);
% derlink = @(mu) exp(mu)./(exp(mu)-1);
% invlink = @(resp) log(1 + exp(resp));
% F = {link, derlink, invlink};

link = @link_test;
derlink = @derlink_test;
invlink = @invlink_test;
F = {link, derlink, invlink};

trainM=[];I_eg=[];ind1st=[];
for ii=1:9
    trainM=[trainM spTrain_byLoc{ii}];
    I_eg=[I_eg Ie_byLoc{ii}];
    ind1st=[ind1st ind1st_allLoc{ii}];
end


%% log-likelihood
gVec=0.02:0.02:0.2;
logL_vec=zeros(length(gVec),1);
for i=1:length(gVec)
    g=gVec(i);
    [expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time
    for iloc=1:8
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


%%
I_eg=repmat(I_e,1,9);
g=gVec(find(logL_vec==max(logL_vec))); %MLE
[expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time
for iloc=1:8
    expg_k_loc(:,iloc)=expg_k(:).*ind_allLoc(:,iloc+1);
end
[betahat_conv_allLoc,~,stats_conv]=glmfit([ind1st(:).*expg_Vreset(:) expg_k(:) expg_k_loc],spTrain_allLoc,'Poisson','link',F);
betahat_conv_allLoc'
[betahat_conv_allLoc(1) E_L-V_th; betahat_conv_allLoc(2) V_reset-E_L;betahat_conv_allLoc(3) k_vec(1);betahat_conv_allLoc(3)+betahat_conv_allLoc(4) k_vec(2);betahat_conv_allLoc(3)+betahat_conv_allLoc(5) k_vec(3);betahat_conv_allLoc(3)+betahat_conv_allLoc(6) k_vec(4);betahat_conv_allLoc(3)+betahat_conv_allLoc(7) k_vec(5);betahat_conv_allLoc(3)+betahat_conv_allLoc(8) k_vec(6);betahat_conv_allLoc(3)+betahat_conv_allLoc(9) k_vec(7);betahat_conv_allLoc(3)+betahat_conv_allLoc(10) k_vec(8);betahat_conv_allLoc(3)+betahat_conv_allLoc(11) k_vec(9);]

%% prediction
locVec=1:9;
figure;suptitle('Prediction (by amplitude)');
simTrainMat=zeros(75,9);
for jj=1:15
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
ntrial=1;
ntime=length(I_eg_fit);
simTrain=zeros(ntime,ntrial);
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
        simTrain(j,tr)=1; %spike
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
        simTrain(j,tr)=0;
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


%% predicted voltage
%%%%for diagnostic purposes, it would be helpful to plot the GLM predicted voltage next to the true spikes,
%%%%especially in cases where the neuron spikes predictably:
%%%%if the predicted voltage consistently gets large in the wrong place, this would indicate that the input current model should be changed

