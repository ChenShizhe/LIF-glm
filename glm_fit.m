%%
%%%%%proper LIF_GLM

I_e=I_e_vect_mat(3:6,2:end)';trainM=spTrain(3:6,:)';

g=0.1;k=1;

%%%%Design matrix
[expg_Vreset,expg_k]=gconv(I_e,trainM,g); %temporally convolve paramters with g upto spike time

%%%%Link function
link = @(mu) log(exp(mu)-1);
derlink = @(mu) exp(mu)./(exp(mu)-1);
invlink = @(resp) log(1 + exp(resp));
F = {link, derlink, invlink};

[betahat_conv,~,stats_conv]=glmfit([expg_Vreset(:) expg_k(:)],trainM(:),'Poisson','link',F);

[betahat_conv(1) -V_th;betahat_conv(2) V_reset;betahat_conv(3) k]







