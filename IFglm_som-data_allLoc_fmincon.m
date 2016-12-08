design = [design1 design2];
num_params=size(design,2);
init_vals=zeros(num_params,1);
x0 = init_vals; x0(1)=-5;x0(2)=-50;

obj_fun = @(x) logL_LIFglmnn( x, design, Y );

% options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','GradObj','on','Diagnostics','on','UseParallel',true); %

options = optimset('Algorithm','interior-point','Display','iter','GradObj','on','Diagnostics','on','UseParallel','always'); %

                            
ub = Inf*ones(size(x0));                            
lb = -Inf*ones(size(x0));
lb(3:num_params) = 0;
ub(3:num_params) = 1e6;

fit_params = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);
