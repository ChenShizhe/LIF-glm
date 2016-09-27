clear all;clc;

%% simulation
figure;suptitle('Simulation (by location)');
for locid=1:9
dt=1; %time step ms
t_end=75;t_StimStart=5;t_StimEnd=15;
V_th=-55; %spike threshold [mV]
E_L=-70; %resting membrane potential [mV]
V_reset=-75; %value to reset voltage to after a spike [mV]
V_spike=20; %value to draw a spike to, when cell spikes
tau=15; %membrane time constant [ms]
R_m=0.375+0.0075*locid; %membrane resistance [MOhm]

stoc_mu=0;stoc_sigma=2;

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
    V_plot_vect(i)=V_vect(i); %if no spike, then just plot the actual voltage V
    V_plot_vect2(i)=V_vect(i);
    I_e_vect=zeros(1,t_StimStart/dt); %portion of I_e_vect from t=0 to t_StimStart
    I_e_vect=[I_e_vect I_Stim*ones(1,1+((t_StimEnd-t_StimStart)/dt))];
    I_e_vect=[I_e_vect zeros(1,(t_end-t_StimEnd)/dt)];
    I_e_vect_mat(:,PlotNum)=I_e_vect;
    
    NumSpikes=0; %holds number of spikes that have occurred
    for t=dt:dt:t_end %loop through values of t in steps of df ms        
        %V_inf = E_L + I_e_vect(i)*R_m; %value that V_vect is exponentially
                                         %decaying towards at this time step
        %V_vect(i+1) = V_inf + (V_vect(i)-V_inf)*exp(-dt/tau);
        
        V_vect(i+1) = V_vect(i) + (E_L-V_vect(i) + I_e_vect(i)*R_m)/tau*dt + sqrt(dt)*normrnd(stoc_mu,stoc_sigma); 
        %Euler's method with stochasticity
        
        
        %if statement below says what to do if voltage crosses threshold
        if (V_vect(i+1)>V_th) %cell spiked
            V_vect(i+1)=V_reset; %set voltage back to V_reset
            V_plot_vect(i+1)=V_spike; %set vector that will be plotted to show a spike here
            V_plot_vect2(i+1) = V_vect(i+1);
            NumSpikes=NumSpikes+1; %add 1 to the total spike count
            spTrain(i,PlotNum)=1;
        else %voltage didn't cross threshold so cell does not spike
            V_plot_vect(i+1)=V_vect(i+1); %plot actual voltage
            V_plot_vect2(i+1)=V_vect(i+1);
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
plot([t_StimStart t_StimStart],[V_reset V_spike],'k:');
plot([t_StimEnd t_StimEnd],[V_reset V_spike],'k:');
plot([0 75],[5.5 5.5],'k-');
plot([0 75],[10.5 10.5],'k-');
plot([0 75],[15.5 15.5],'k-');
hold off
if (locid==1 | locid==4 | locid==7)
    ylabel('25mV 50mV 100mV');
end
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
    set(gca,'FontSize',12);
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


%% glm fit

spTrain_allLoc=[];
Ie_allLoc=[];
ind_allLoc=zeros(size(spTrain,1)*size(spTrain,2),9);
for ii=1:9
    spTrain=spTrain_byLoc{ii};
    spTrain_allLoc=[spTrain_allLoc;spTrain(:)];
    I_e=Ie_byLoc{ii};
    Ie_allLoc=[Ie_allLoc;I_e(:)];
    nii=size(spTrain,1)*size(spTrain,2);
    ind_allLoc((ii-1)*nii+1:ii*nii,ii)=ii-1;
end


%% logLikelhood
link = @(mu) log(exp(mu)-1);
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

trainM=[];I_eg=[];
for ii=1:9
    trainM=[trainM spTrain_byLoc{ii}];
    I_eg=[I_eg Ie_byLoc{ii}];
end

gVec=0.02:0.02:1;
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


%%
I_eg=repmat(I_e,1,9);
g=gVec(find(logL_vec==max(logL_vec)));
[expg_Vreset,expg_EL,expg_k]=gconv(I_eg,trainM,g); %temporally convolve paramters with g upto spike time

%%%%all loc with an indicator function
[betahat_conv_allLoc,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:) ind_allLoc(:,2:9)],spTrain_allLoc,'Poisson','link',F);
betahat_conv_allLoc'


%% prediction
locVec=1:9;
figure;suptitle('Prediction (by amplitude)');
simTrainMat=zeros(75,9);clear betahat_conv
for jj=1:15
for ii=1:9
    betahat_conv=zeros(1,4);
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


