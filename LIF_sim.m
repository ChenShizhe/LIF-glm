%%
%%%%simulation
clear all;
close all;

%%%DEFINE PARAMETERS
dt=1; %time step ms
t_end=500; %total run time ms
t_StimStart=100; %time to start injecting current
t_StimEnd=400; %time to end injecting current
V_th=-55; %spike threshold
V_reset=-70; %value to reset voltage to after a spike
V_spike=20; %value to draw a spike to, when cell spikes
k=1;g=0.1;

%%%DEFINE INITIAL VALUES AND VECTORS TO HOLD RESULTS
t_vect=0:dt:t_end; 
V_vect=zeros(1,length(t_vect));
V_plot_vect=zeros(1,length(t_vect));

%INTEGRATE THE EQUATION dV/dt = -g*(V-V_reset)+k*I_e
PlotNum=0;
I_Stim_vect=1.43:0.04:1.63; %magnitudes of pulse of injected current [nA]
spTrain=zeros(t_end,length(I_Stim_vect));

for I_Stim=I_Stim_vect; %loop over different I_Stim values
    PlotNum=PlotNum+1;
    i=1; %index denoting which element of V is being assigned
    V_vect(i)=V_reset; %first element of V, i.e. value of V at t=0
    V_plot_vect(i)=V_vect(i); %if no spike, then just plot the actual voltage V
    I_e_vect=zeros(1,t_StimStart/dt); %portion of I_e_vect from t=0 to t_StimStart
    I_e_vect=[I_e_vect I_Stim*ones(1,1+((t_StimEnd-t_StimStart)/dt))];
    I_e_vect=[I_e_vect zeros(1,(t_end-t_StimEnd)/dt)];
    I_e_vect_mat(:,PlotNum)=I_e_vect;
    
    NumSpikes=0; %holds number of spikes that have occurred
    for t=dt:dt:t_end %loop through values of t in steps of df ms        
        dV=(-g*(V_vect(i)-V_reset)+k*I_e_vect(i))*dt; %Euler's method
        V_vect(i+1)=V_vect(i)+dV;
        
        %if statement below says what to do if voltage crosses threshold
        if (V_vect(i+1)>V_th) %cell spiked
            V_vect(i+1)=V_reset; %set voltage back to V_reset
            V_plot_vect(i+1)=V_spike; %set vector that will be plotted to show a spike here
            NumSpikes=NumSpikes+1; %add 1 to the total spike count
            spTrain(t,PlotNum)=1;
        else %voltage didn't cross threshold so cell does not spike
            V_plot_vect(i+1)=V_vect(i+1); %plot actual voltage
        end
        i=i+1; %add 1 to index,corresponding to moving forward 1 time step
    end    
    
    
    %MAKE PLOTS
    figure(1)
    subplot(length(I_Stim_vect),1,PlotNum)
    plot(t_vect,V_plot_vect);
    if (PlotNum==1)
        title('Voltage vs. time');
    end
    if (PlotNum==length(I_Stim_vect))
        xlabel('Time in ms');
    end
    ylabel('Voltage in mV');
end
