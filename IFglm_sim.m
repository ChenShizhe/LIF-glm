clear all
%%%DEFINE PARAMETERS
dt=1; %time step ms
t_end=75; %total run time ms
t_StimStart=6; %time to start injecting current
t_StimEnd=15; %time to end injecting current
V_th=-20; %spike threshold [mV]
E_L=-26; %resting membrane potential [mV]
V_reset=-75; %value to reset voltage to after a spike [mV]
V_spike=20; %value to draw a spike to, when cell spikes
tau=5; %membrane time constant [ms]
R_m=0.178; %membrane resistance [MOhm]

%%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect=0:dt:t_end; 
V_vect=zeros(1,length(t_vect));
V_plot_vect=zeros(1,length(t_vect));

%INTEGRATE THE EQUATION tau*dV/dt = -V + E_L + I_e*R_m, I_e with Gaussian noise
PlotNum=0;
I_Stim_vect=20:20:100;
%I_Stim_vect=reshape(repmat(20:20:100,5,1),1,25); %magnitudes of pulse of injected current [nA]
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
        %V_inf = E_L + I_e_vect(i)*R_m;
        %V_vect(i+1) = V_inf + (V_vect(i)-V_inf)*exp(-dt/tau);
        
        V_vect(i+1) = V_vect(i) + (E_L-V_vect(i) + I_e_vect(i)*R_m)/tau; %Euler's method
        
        
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
    
    
    %MAKE PLOTS
    figure(1)
    subplot(length(I_Stim_vect),1,PlotNum)
    plot(t_vect,V_plot_vect);
    xlim([0 t_end]);
    ylim([V_reset V_spike]);
    hold on
    plot([t_StimStart t_StimStart],[V_reset V_spike],'k:');
    plot([t_StimEnd t_StimEnd],[V_reset V_spike],'k-');
    hold off
    if (PlotNum==1)
        title(['g=' num2str(1/tau) 'k=' num2str(R_m/tau) 'E\_L=' num2str(E_L) ', V\_res=' num2str(V_reset) ', V\_th=' num2str(V_th)]);
    end
    if (PlotNum==length(I_Stim_vect))
        xlabel('Time (ms)');
    end
    ylabel('Voltage (mV)');
end


%%
k=R_m/tau;g=1/tau;

I_e=I_e_vect_mat(2:end,:);trainM=spTrain;

[expg_Vreset,expg_EL,expg_k]=gconv(I_e,trainM,g); %temporally convolve paramters with g upto spike time

figure;
plot(t_vect(2:end),V_reset.*expg_Vreset(:,end) +E_L.*expg_EL(:,end) + k.*expg_k(:,end),'r');
set(gca,'FontSize',16);
xlabel('Time (ms)');ylabel('Voltage (mV)');
hold on
plot(t_vect(2:end),V_plot_vect2(2:end),'b');
plot([0 t_end],[V_th V_th],'k');
scatter(find(trainM(:,end)),V_vect(find(trainM(:,end))+1))
hold off
xlim([0 t_end]);
legend('Fitted V(t)','True V(t)');

%%
link = @(mu) log(exp(mu)-1);
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:)],trainM(:),'Poisson','link',F);

[betahat_conv(1) E_L-V_th;betahat_conv(2) V_reset-E_L;betahat_conv(3) k]

%% logLikelihood
I_e=I_e_vect_mat(2:end,:);trainM=spTrain;
gVec=0.02:0.04:0.3;logL_vec=zeros(length(gVec),1);
for i=1:length(gVec)
    g=gVec(i);

clear expg_Vreset expg_EL expg_k
[expg_Vreset,expg_EL,expg_k]=gconv(I_e,trainM,g);

[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:)],trainM(:),'Poisson','link',F);
lambdahat=glmval(betahat_conv,[expg_Vreset(:) expg_k(:)],F);
logL=sum(log(poisspdf(trainM(:),lambdahat)));

logL_vec(i)=logL;
end

figure;plot(logL_vec);set(gca,'XTick',1:length(gVec),'XTickLabel',gVec,'FontSize',16);
xlabel('g');ylabel('log-likelihood');
[find(logL_vec==max(logL_vec)) gVec(find(logL_vec==max(logL_vec)))]



%% prediction
%nI=size(I_e_vect_mat,2);
g=1/tau;
[expg_Vreset,expg_EL,expg_k]=gconv(I_e,trainM,g); %temporally convolve paramters with g upto spike time

[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:)],trainM(:),'Poisson','link',F);
betahat_conv'

I_eg_fit=I_e(:,end);
ntrial=20;
ntime=length(I_eg_fit);
simTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);
spiketime=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);
fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);
    j=1;lastSpT=0;
    fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
    fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
    step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr))+1);
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
            step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr))+1);
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
            step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr))+1);
            lambdaMat(j,tr)=step(j);
        end
    end
end
clear tao step
end

figure;
subplot(5,1,1);
stem(t_vect(2:end),spTrain(:,end));
set(gca,'FontSize',12);title('Simulation (Stim 100 mV)');
%imagesc(t_vect(2:end),1:5,spTrain(:,end-4:end)');set(gca,'FontSize',12);
%colormap(flipud(gray));
%[r1,n1]=max(spTrain~=0,[],1);
box off;axis off;
%title({'Simulation',num2str([round(mean(sum(spTrain,1))*100)/100 round(std(sum(spTrain,1))*100)/100 mean(n1)])});
subplot(5,1,[2:5]);
imagesc(simTrain');set(gca,'FontSize',12);
colormap(flipud(gray));
xlabel('Time (ms)');ylabel('Trial');
[r2,n2]=max(simTrain~=0,[],1);
box off;
title({'Prediction',num2str([round(mean(sum(simTrain,1))*100)/100 round(std(sum(simTrain,1))*100)/100 mean(n2)])});







% % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% with an indicator function, 1 = not yet 1st spike %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % %



%%
ind1stM=zeros(size(trainM));
for i=1:size(trainM,2)
    sp1st=find(trainM(:,i),1);
    ind1stM(1:sp1st,i)=1;
end
[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:) ind1stM(:)],trainM(:),'Poisson','link',F);
[betahat_conv(1) E_L-V_th;betahat_conv(2) V_reset-E_L;betahat_conv(3) k]

lambdahat=glmval(betahat_conv,[expg_Vreset(:) expg_k(:)  ind1stM(:)],F);
logL=sum(log(poisspdf(trainM(:),lambdahat)));

%% logLikelihood
I_e=I_e_vect_mat(2:end,:);trainM=spTrain;
ind1stM=zeros(size(trainM));
for i=1:size(trainM,2)
    sp1st=find(trainM(:,i),1);
    ind1stM(1:sp1st,i)=1;
end

gVec=0.02:0.04:0.4;logL_vec=zeros(length(gVec),1);
for i=1:length(gVec)
    g=gVec(i);

clear expg_Vreset expg_EL expg_k
[expg_Vreset,expg_EL,expg_k]=gconv(I_e,trainM,g);


[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:) ind1stM(:)],trainM(:),'Poisson','link',F);
lambdahat=glmval(betahat_conv,[expg_Vreset(:) expg_k(:)  ind1stM(:)],F);

%[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:)],trainM(:),'Poisson','link',F);
%lambdahat=glmval(betahat_conv,[expg_Vreset(:) expg_k(:)],F);

logL=sum(log(poisspdf(trainM(:),lambdahat)));

logL_vec(i)=logL;
end

figure;plot(logL_vec);set(gca,'XTick',1:length(gVec),'XTickLabel',gVec,'FontSize',16);
xlabel('g');ylabel('log-likelihood');
[find(logL_vec==max(logL_vec)) gVec(find(logL_vec==max(logL_vec)))]



%% prediction
%nI=size(I_e_vect_mat,2);
g=1/tau;
[expg_Vreset,expg_EL,expg_k]=gconv(I_e,trainM,g); %temporally convolve paramters with g upto spike time

[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:) ind1stM(:)],trainM(:),'Poisson','link',F);
betahat_conv'

I_eg_fit=I_e(:,end);
ntrial=20;
ntime=length(I_eg_fit);
simTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);
spiketime=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);
fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);
    j=1;lastSpT=0;
    fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
    fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*I_eg_fit(lastSpT+1:j);
    step(j)=log(exp(betahat_conv(1)+betahat_conv(4)*min(1,lastSpT+1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr))+1);
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
            step(j)=log(exp(betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr))+1);
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
            step(j)=log(exp(betahat_conv(1)+betahat_conv(4)*min(1,lastSpT+1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr))+1);
            lambdaMat(j,tr)=step(j);
        end
    end
end
clear tao step
end

figure;
subplot(5,1,1);
stem(t_vect(2:end),spTrain(:,end));
set(gca,'FontSize',12);title('Simulation (Stim 100 mV)');
%imagesc(t_vect(2:end),1:5,spTrain(:,end-4:end)');set(gca,'FontSize',12);
%colormap(flipud(gray));
%[r1,n1]=max(spTrain~=0,[],1);
box off;axis off;
%title({'Simulation',num2str([round(mean(sum(spTrain,1))*100)/100 round(std(sum(spTrain,1))*100)/100 mean(n1)])});
subplot(5,1,[2:5]);
imagesc(simTrain');set(gca,'FontSize',12);
colormap(flipud(gray));
xlabel('Time (ms)');ylabel('Trial');
[r2,n2]=max(simTrain~=0,[],1);
box off;
title({'Prediction',num2str([round(mean(sum(simTrain,1))*100)/100 round(std(sum(simTrain,1))*100)/100 mean(n2)])});



