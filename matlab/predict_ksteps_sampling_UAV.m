function [pred_X]=predict_ksteps_sampling_UAV(N,dOmega,w,para,x0,dynamic,pre_kstep)
   %N: the number of testing samples of trajectory  
    %prior=para.prior;
    sam_omega=w;
    pred_X=zeros(size(x0,1),N,pre_kstep);
    cur_X=repmat(x0,1,N);
    mu=zeros(size(x0,1),1);
    for k=1:pre_kstep
        for i=1:N

            dOmega_s=dOmega(1)+(dOmega(2)-dOmega(1)).*rand(); 
            while(sam_omega+dOmega_s>para.max ||sam_omega+dOmega_s<para.min)
                dOmega_s=dOmega(1)+(dOmega(2)-dOmega(1)).*rand(); 
            end
            sam_omega = sam_omega+dOmega_s;
            
            pred_X(:,i,k) = dynamic.fun(sam_omega,cur_X(:,i),dynamic.h);%+ mvnrnd(mu,CT.cov,1)';
            cur_X(:,i) = pred_X(:,i,k);
        end
        %prior = (1-para.epsilon)*prior +para.epsilon*para.prior_initial;
    end
end