function [post]=computePosterior_Dubinsmodel(para,x,next_x,dynamic) 
 %%compute posterior p(omega(t)|x(t+1))
        for i=1:para.num_p
            
            xmean=dynamic.fun(para.theta(i,:),x,dynamic.h);
           
            p=mvnpdf(xmean(:,1),next_x,dynamic.cov);

            post(i)=p*para.prior(i);
            
        end
        normC=sum(post);
        post=post/normC;
end