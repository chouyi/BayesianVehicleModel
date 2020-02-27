function [post]=computePosterior_CTmodel(para,x,next_x,CT) 
 %%compute posterior p(omega(t)|x(t+1))
        for i=1:para.num_p
            xmean=CT.A(:,:,i) * x;
           
            p=mvnpdf(xmean(:,1),next_x,CT.cov);

            post(i)=p*para.prior(i);
            
        end
        normC=sum(post);
        post=post/normC;
end