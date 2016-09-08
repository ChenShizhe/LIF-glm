load('data\spikes_fix.mat');

%%
ce=1; %choose a cell

for bl=1:15
    spBlockAll1ms{bl}=spikes_downres{ce}(1+(bl-1)*121:121+(bl-1)*121,:);
end

for bl=15
    spBlock=spBlockAll1ms{bl};
spO=zeros(121,1);
for loc=1:121
    if isempty(find(spBlock(loc,:)))==0
        spO(loc)=find(spBlock(loc,:),1);
    end
end
[s1 s2]=sort(spO);
locAll{bl}=s2(find(s1)); %locations that spike

figure;
imagesc(spBlock(s2,:));%title(['Block ',int2str(bl)]);
set(gca,'YTick',1:1:121,'YTickLabel',s2,'XTick',0:5:75,'XTickLabel',0:5:75);
colormap(flipud(gray));
xlabel('Time (ms)','FontSize',20);ylabel('Location','FontSize',20);
box off;
hold on
plot([5 5],[0 122],'k:');
plot([15 15],[0 122],'k:');
hold off
xlim([0 75]);ylim([75 121]);
end

%%
load('data\current_template.mat');

loc=61;

sp_eg_100mV=[spBlockAll1ms{11}(loc,:)' spBlockAll1ms{12}(loc,:)' spBlockAll1ms{13}(loc,:)' spBlockAll1ms{14}(loc,:)' spBlockAll1ms{14}(loc,:)'];
sp_eg_50mV=[spBlockAll1ms{6}(loc,:)' spBlockAll1ms{7}(loc,:)' spBlockAll1ms{8}(loc,:)' spBlockAll1ms{9}(loc,:)' spBlockAll1ms{10}(loc,:)'];
sp_eg_25mV=[spBlockAll1ms{1}(loc,:)' spBlockAll1ms{2}(loc,:)' spBlockAll1ms{3}(loc,:)' spBlockAll1ms{4}(loc,:)' spBlockAll1ms{5}(loc,:)'];

I_eg_100mV=norm_average_current(1:20:20*75)*100;
I_eg_50mV=norm_average_current(1:20:20*75)*50;
I_eg_25mV=norm_average_current(1:20:20*75)*25;

I_e=[repmat(I_eg_100mV',1,5) repmat(I_eg_50mV',1,5) repmat(I_eg_25mV',1,5)];
spTrain=[sp_eg_100mV sp_eg_50mV sp_eg_25mV];


%%
%%%%dV/dt = -g(V_t-E_L) + k*x(t)
%%%%lambda_t = f(V_t)
%%%%V_t  = (E_leak - V_thres)
%%%%      + (V_reset - E_leak)*exp(-gt)
%%%%      + k*int I(s)*exp(-g(t-s))ds

g=0.1;

clear expg_Vreset expg_k
[expg_Vreset,expg_k]=gconv(I_e,spTrain,g); %temporally convolve paramters with g upto spike time

link = @(mu) log(exp(mu)-1);
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:)],spTrain(:),'Poisson','link',F);
betahat_conv'

figure;
plot(betahat_conv(1)+betahat_conv(2).*expg_Vreset(:)+betahat_conv(3).*expg_k(:),'r');
set(gca,'FontSize',16);
xlabel('Time (ms)');ylabel('Voltage (mV)');
hold on
plot(t_vect(2:end),V_plot_vect2(2:end),'b');
plot([0 t_end],[V_th V_th],'k');
scatter(find(trainM(:,end)),V_vect(find(trainM(:,end))+1))
% plot(t_vect(2:end),lambda(2:end)*50 - 65)
hold off
xlim([0 t_end]);
legend('Fitted V(t)','True V(t)');




%% prediction
I_eg_fit=I_eg_100mV;
betahat=betahat_conv;
ntrial=10;
ntime=75;
simTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);
spiketime=zeros(ntime,ntrial);
spikebin=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);
fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);
    j=1;
    t_elapse=j;
    fit_expg_Vreset(j,tr)=exp(-g*t_elapse);
    fit_expg_k(j,tr)=exp(-g*t_elapse)*I_eg_fit(j);
    a=[fit_expg_Vreset(j,tr) sum(fit_expg_k(1:j,tr))]; 
    step(j)=betahat(1)+dot(betahat(2:end),a);
    lambdaMat(j,tr)=step(j);
    spikebin(:,tr)=0;
while (j<=ntime)
    if sum(step)>tao
        sim_Zi(j,tr)=tao;
        simTrain(j,tr)=1; %spike
        spikebin(:,tr)=j;
        t_elapse=0;
        if j==ntime
            break
        else
            j=j+1;
            t_elapse=t_elapse+1;
            tao=exprnd(1);
            step(1:j-1)=0;
            fit_expg_k(1:j-1,tr)=0;
            fit_expg_k(j,tr)=exp(-g*t_elapse)*I_eg_fit(j); %restart rolling sum for parameter k
            fit_expg_Vreset(j,tr)=exp(-g*t_elapse);
            a=[fit_expg_Vreset(j,tr) sum(fit_expg_k(1:j,tr))]; 
            step(j)=betahat(1)+dot(betahat(2:end),a);
            lambdaMat(j,tr)=step(j);
        end
    else
        simTrain(j,tr)=0;
        if j==ntime
            break
        else
            j=j+1;
            t_elapse=t_elapse+1;
            fit_expg_k(j,tr)=exp(-g*t_elapse)*I_eg_fit(j); %keep rolling sum for parameter k
            fit_expg_Vreset(j,tr)=exp(-g*t_elapse);
            a=[fit_expg_Vreset(j,tr) sum(fit_expg_k(1:j,tr))]; 
            step(j)=betahat(1)+dot(betahat(2:end),a);
            lambdaMat(j,tr)=step(j);
        end
    end
end
clear step tao
end

figure;
imagesc(simTrain');
colormap(flipud(gray));
xlabel('Time (ms)');ylabel('Trial (100mV)');
box off;



%%





%%


