%% Simulation

%%%DEFINE PARAMETERS
dt=1; %time step ms
t_end=75; %total run time ms
t_StimStart=6; %time to start injecting current
t_StimEnd=15; %time to end injecting current
V_th=-55; %spike threshold [mV]
E_L=-60; %resting membrane potential [mV]
V_reset=-60; %value to reset voltage to after a spike [mV]
V_spike=20; %value to draw a spike to, when cell spikes
tau=15; %membrane time constant [ms]
R_m=10; %membrane resistance [MOhm]

%%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect=0:dt:t_end; 
V_vect=zeros(1,length(t_vect));
V_plot_vect=zeros(1,length(t_vect));
V_plot_vect2=zeros(1,length(t_vect));


%%
%INTEGRATE THE EQUATION tau*dV/dt = -V + E_L + I_e*R_m
PlotNum=0;
I_Stim_vect=.5:0.1:1.0; %magnitudes of pulse of injected current [nA]
spTrain=zeros(t_end,length(I_Stim_vect));
clear I_e_vect_mat

for I_Stim=I_Stim_vect; %loop over different I_Stim values
    
    PlotNum = PlotNum+1;
    i=1; %index denoting which element of V is being assigned
        
    V_vect(i)=E_L; %first element of V, i.e. value of V at t=0
    V_plot_vect(i)=V_vect(i); %if no spike, then just plot the actual voltage V
    V_plot_vect2(i)=V_vect(i);
    I_e_vect=zeros(1,t_StimStart/dt); %portion of I_e_vect from t=0 to t_StimStart
    I_e_vect=[I_e_vect I_Stim*ones(1,1+((t_StimEnd-t_StimStart)/dt))];
    I_e_vect=[I_e_vect zeros(1,(t_end-t_StimEnd)/dt)];
%     I_e_vect=awgn(I_e_vect,10,'measured');
    I_e_vect_mat(:,PlotNum)=I_e_vect;
    
    NumSpikes=0; %holds number of spikes that have occurred
    
    tao=exprnd(1);
    lambda(i)=log(exp(V_vect(i)-V_th)+1);
    last_spike = 1;
    
    for t=dt:dt:t_end %loop through values of t in steps of df ms        
        %V_inf = E_L + I_e_vect(i)*R_m;
        %V_vect(i+1) = V_inf + (V_vect(i)-V_inf)*exp(-dt/tau);
        
        V_vect(i+1) = V_vect(i) + (E_L-V_vect(i) + I_e_vect(i)*R_m*5)/tau; %Euler's method
        lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
        
        %if statement below says what to do if voltage crosses threshold
        if sum(lambda(last_spike+1:i+1))>tao %cell spiked
            V_plot_vect2(i+1) = V_vect(i+1);
            V_vect(i+1)=V_reset; %set voltage back to V_reset
            V_plot_vect(i+1)=V_spike; %set vector that will be plotted to show a spike here
            
            NumSpikes=NumSpikes+1; %add 1 to the total spike count
            spTrain(i,PlotNum)=1;
            if last_spike == i+1
                disp('two spikes in a row!!')
                return
            end
            last_spike = i+1;
            tao=exprnd(1);
%             lambda(1:i)=0;
%             lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
        else %voltage didn't cross threshold so cell does not spike
            V_plot_vect(i+1)=V_vect(i+1); %plot actual voltage
            V_plot_vect2(i+1) = V_vect(i+1);
%             lambda(i+1)=log(exp(V_vect(i+1)-V_th)+1);
        end
        i=i+1; %add 1 to index,corresponding to moving forward 1 time step
    end    
    
    
    %MAKE PLOTS
    figure(2)
    subplot(length(I_Stim_vect),1,PlotNum)
    plot(t_vect,V_plot_vect);
%     hold on
%     plot(t_vect,lambda*100)
    if (PlotNum==1)
        title('Voltage vs. time');
    end
    if (PlotNum==length(I_Stim_vect))
        xlabel('Time (ms)');
    end
    ylabel('Voltage (mV)');
    ylim([-100 100])
    
    if I_Stim ~= I_Stim_vect(end)
        clear lambda
        V_vect=zeros(1,length(t_vect));
        V_plot_vect=zeros(1,length(t_vect));
    end
end


%%
k=R_m/tau;g=1/tau;

I_e=I_e_vect_mat(1:end-1,:);trainM=spTrain;

[expg_Vreset,expg_EL,expg_k]=gconv(I_e,trainM,g); %temporally convolve paramters with g upto spike time

figure;
plot(t_vect(2:end),V_reset.*expg_Vreset(:,end) +E_L.*expg_EL(:,end) + k.*expg_k(:,end),'r');
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

%%
link = @(mu) log(exp(mu)-1);  %link = @(mu) mu + log(1-exp(-mu));
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:)],trainM(:),'Poisson','link',F);

[betahat_conv(1) E_L-V_th;betahat_conv(2) V_reset-E_L;betahat_conv(3) k]

