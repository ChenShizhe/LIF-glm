%%
%%%%%proper LIF_GLM

I_e=I_e_vect_mat(2:end,:);trainM=spTrain;

[expg_Vreset,expg_EL,expg_k]=gconv(I_e,trainM,g); %temporally convolve paramters with g upto spike time

figure;
plot(t_vect(2:end),V_reset.*expg_Vreset(:,end) +E_L.*expg_EL(:,end) + k.*expg_k(:,end),'r.');
hold on
plot(t_vect(2:end),V_vect(2:end),'b.');
plot([0 t_end],[-55 -55],'k');
hold off
xlim([0 t_end]);
legend('Fitted V(t)','True V(t)');


%%%%collinear
figure;
plot(expg_Vreset(:,end), expg_k(:,end),'.');
