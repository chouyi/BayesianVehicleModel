function [pred_X]=predict_ksteps_sampling_Dubins(N,w,para,x0,dynamic,pre_kstep)
   %N: the number of testing samples of trajectory  
    %prior=para.prior;
    theta=w;
    
    pred_X=zeros(size(x0,1),N,pre_kstep);
    cur_X=repmat(x0,1,N);
    %mu=zeros(size(x0,1),1);
    for k=1:pre_kstep
        for i=1:N
%             max_theta=theta+para.dpmax;
%             min_theta=theta+para.dpmin;
%             if(max_theta(1)>para.max(1))
%                 max_theta(1)=para.max(1);
%             end
%             if(max_theta(2)>para.max(2))
%                 max_theta(2)=para.max(2);
%             end
%             if(min_theta(1)<para.min(1))
%                 min_theta(1)=para.min(1);
%             end
%             if(min_theta(2)<para.min(2))
%                 min_theta(2)=para.min(2);
%             end
%             theta=min_theta+(max_theta-min_theta).*rand(1,para.np);


            dTheta=para.dpmin+(para.dpmax-para.dpmin).*rand(1,para.np);

            while(any(theta+dTheta>para.max) ||any(theta+dTheta<para.min))
                dTheta=para.dpmin+(para.dpmax-para.dpmin).*rand(1,para.np); 
            end
            theta = theta+dTheta;
            pred_X(:,i,k) = dubinsDynamic(theta,cur_X(:,i),dynamic.h);%+ mvnrnd(mu,model_cov,1)';
            %pred_X(:,i,k) = dynamic.fun(sam_omega,cur_X(:,i),dynamic.h);%+ mvnrnd(mu,CT.cov,1)';
            cur_X(:,i) = pred_X(:,i,k);
        end
        %prior = (1-para.epsilon)*prior +para.epsilon*para.prior_initial;
    end
end