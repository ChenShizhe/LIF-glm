%%
%%%%%proper LIF_GLM

%%
I_e=I_e_vect_mat(2:end,:);trainM=spTrain;

[expg_Vreset,expg_k]=gconv(I_e,trainM,g); %temporally convolve paramters with g upto spike time

figure;
subplot(511);
plot(V_reset.*expg_Vreset(:,end));
%ylim([0 1]);
set(gca,'FontSize',12);
ylabel('Convoluted g');
subplot(512);
plot(k.*expg_k(:,end));
set(gca,'FontSize',12);
ylabel('Convoluted I');
subplot(513);
plot(V_reset.*expg_Vreset(:,end) + k.*expg_k(:,end));
ylim([-70 -50]);
set(gca,'FontSize',12);
ylabel('Fitted V(t)');
subplot(514);
plot(t_vect,V_vect);
set(gca,'FontSize',12);
ylabel('True V(t)');
subplot(515);
plot(trainM(:,end));
set(gca,'FontSize',12);
ylabel('True spikes');


%%
%%%%glm fit
link = @(mu) log(exp(mu)-1);
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:)],trainM(:),'Poisson','link',F);

[betahat_conv(1) -V_th;betahat_conv(2) V_reset;betahat_conv(3) k]






