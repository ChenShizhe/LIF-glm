%% proper LIF_GLM
k=R_m/tau;g=1/tau;

I_e=I_e_vect_mat(1:end-1,:);trainM=spTrain;

[expg_Vreset,expg_EL,expg_k]=gconv(I_e,trainM,g); %temporally convolve paramters with g upto spike time

figure;
plot(t_vect(2:end),V_reset.*expg_Vreset(:,end) +E_L.*expg_EL(:,end) + k.*expg_k(:,end),'r.');
set(gca,'FontSize',16);
xlabel('Time (ms)');ylabel('Voltage (mV)');
hold on
plot(t_vect(2:end),V_vect(2:end),'b.');
plot([0 t_end],[V_th V_th],'k');
hold off
xlim([0 t_end]);
legend('Fitted V(t)','True V(t)');

%%%%check for collinearity in the design matrix
figure;
plot(expg_Vreset(:,end), expg_k(:,end),'.');
set(gca,'FontSize',16);
xlabel('Design column for V_{reset}');ylabel('Design column for k');


%%%%approximate the inverse link function max(0, V(t)-V_thres) with a soft rectifying function f(lambda)=log(1+exp(lambda))
link = @(mu) log(exp(mu)-1);  %link = @(mu) mu + log(1-exp(-mu));
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:)],trainM(:),'Poisson','link',F);

[betahat_conv(1) E_L-V_th;betahat_conv(2) V_reset-E_L;betahat_conv(3) k]
