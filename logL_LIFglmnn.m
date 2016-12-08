function [ log_likelihood, gradient_ll ] = logL_LIFglmnn( x, design, Y )

num_params=size(design,2);
num_T=size(design,1);

%% log-likelihood
step=zeros(num_T,1);
lambda=zeros(num_T,1);
for i=1:num_T
    step(i) = sum(x.*design(i,:));
%     if step(i)>0
%         lambda(i)=step(i)+1;
%     else
%         lambda(i)=exp(step(i));
%     end
    if step(i)>0
        lambda(i)=1;
    else
        lambda(i)=0;
    end
end

log_likelihood = sum(sum(-step.*Y + lambda));

%% gradient

grad_ll = zeros(num_params,1);
x = ones(size(Y));

grad_ll(1) = sum(sum(-design(:,1).*Y + lambda.*design(:,1)));
grad_ll(2) = sum(sum(-design(:,2).*Y + lambda.*design(:,2)));
for i=3:num_params
grad_ll(i)=sum(sum(-design(:,i).*Y + lambda.*design(:,i)));
end

end
