%%compute the prob of collision by sampling
%% input:sampling states,obstacle
%% output: probability
function [prob] = Prob_collision(X,obs)
    Ns=size(X,2);
    sum1=0;
    for i=1:Ns
        if(X(1,i)>=obs(1,1) && X(1,i)<=obs(1,2)&&X(3,i)>=obs(1,3) && X(3,i)<=obs(1,4))
            sum1=sum1+1;
        end
    end
    prob=sum1/Ns;

end