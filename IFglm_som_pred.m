function [ pred_meanSp, pred_stdSp, pred_latency1st ] = IFglm_som_pred( nloc_sp, betahat_conv_allLoc, ind_e, e_max, all_data )
%LIF-glmnet prediction for som data, matlab portion
%%% input:
%%%1. nloc_sp: number of location;
%%%2. betahat_conv_allLoc: parameter estimates from glmnet in R;
%%%3. ind_e: spike train ID from Excel;
%%%4. e_max: if 1, use input current at max loc only; if 0, use all.
%%%5. all_data: som data set.

gain_vec=zeros(nloc_sp,1);
gain_vec(loc_min)=betahat_conv_allLoc(3);
loc_mmVec=[1:loc_min-1 loc_min+1:nloc_sp];
for i=1:nloc_sp-1
    gain_vec(loc_mmVec(i))=betahat_conv_allLoc(3)+betahat_conv_allLoc(3+i);
end

%%

locVec=1:25;Ie_byLoc=[];
for loc=locVec
    Ie=repmat(all_data{ind_e}{loc}(:,1:20:2000)',1,9)*-1;
    Ie_byLoc{loc}=Ie;
end

ntrial=100; perloc=9;
pred_meanSp=zeros(1,nloc);
pred_stdSp=zeros(1,nloc);
pred_latency1st=zeros(1,nloc); %matrix of mean latency for 1st spike: amp*loc
for ii=1:nloc_sp
    loc=locVec_sp(ii);
    
    if e_max==0
        loc_e=loc;
    elseif e_max==1
        loc_e=12;
    end
    
    Ie_fit=Ie_byLoc{loc_e}(:,1);
    
for jj=1
    betahat_conv=zeros(1,3);
    betahat_conv(1:2)=betahat_conv_allLoc(1:2);
    betahat_conv(3)=gain_vec(ii);
    
ntime=length(Ie_fit);
predTrain=zeros(ntime,ntrial);
lambdaMat=zeros(ntime,ntrial);spiketime=zeros(ntime,ntrial);
fit_expg_Vreset=zeros(ntime,ntrial);fit_expg_k=zeros(ntime,ntrial);
for tr=1:ntrial
    tao=exprnd(1);step=zeros(ntime,1);
    j=1;lastSpT=0;numSp=0;
    fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
    fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*Ie_fit(lastSpT+1:j);
    
    %%%%link function f(x) = 1+x, x>0; = exp(x), x<0.
    
    step0=betahat_conv(1)+betahat_conv(3)*fit_expg_k(j,tr);
    if step0>0
        step(j)=step0+1;
    else
        step(j)=exp(step0);
    end
    
    lambdaMat(j,tr)=step(j);
    
while (j<=ntime)
    if sum(step)>tao
        step(1:j)=0;
        predTrain(j,tr)=1; %spike
        tao=exprnd(1);
        lastSpT=j;numSp=numSp+1;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*Ie_fit(lastSpT+1:j);
            
            step0=betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr);
            if step0>0
                step(j)=step0+1;
            else
                step(j)=exp(step0);
            end
            
            lambdaMat(j,tr)=step(j);
        end
    else
        predTrain(j,tr)=0;
        if j==ntime
            break
        else
            j=j+1;
            fit_expg_Vreset(j,tr)=exp(-g*(j-lastSpT));
            fit_expg_k(j,tr)=exp(-g.*(j-[lastSpT+1:j]))*Ie_fit(lastSpT+1:j);
            
            if numSp==0
                step0=betahat_conv(1)+betahat_conv(3)*fit_expg_k(j,tr);
                if step0>0
                    step(j)=step0+1;
                else
                    step(j)=exp(step0);
                end
            else
                step0=betahat_conv(1)+betahat_conv(2)*fit_expg_Vreset(j,tr)+betahat_conv(3)*fit_expg_k(j,tr);
                if step0>0
                    step(j)=step0+1;
                else
                    step(j)=exp(step0);
                end
            end
            
            lambdaMat(j,tr)=step(j);
        end
    end
end
end

[r2,n2]=max(predTrain~=0,[],1);n2(r2==0)=ntime;pred_latency1st(jj,loc)=mean(n2);
pred_meanSp(jj,loc)=round(mean(sum(predTrain,1))*100)/100;
pred_stdSp(jj,loc)=round(std(sum(predTrain,1))*100)/100;
end
predTrain_allLoc{loc}=predTrain;
end
pred_latency1st(pred_latency1st==0)=ntime;

end

