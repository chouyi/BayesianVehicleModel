% %%compute the prob of collision by sampling
% %% input:sampling states,obstacle
% %% output: probability
% function [prob] = Prob_collision(X,Y,obs)
%     Ns=size(X,2);
%     sum1=0;
%     for i=1:Ns
%         if(X(1,i)>=obs(1,1) && X(1,i)<=obs(1,2)&&Y(1,i)>=obs(1,3) && Y(1,i)<=obs(1,4))
%             sum1=sum1+1;
%         end
%     end
%     prob=sum1/Ns;
% 
% end

function [prob] = Prob_collision(x, y,obs_zone)
    Ns=size(x,2);
    sum1=0;
    for i=1:Ns
        collision = inpolygon(x(1, i), y(1,i), obs_zone.Vertices(:,1),obs_zone.Vertices(:,2));
        if collision == 1
            sum1=sum1 + 1;
        end
    end
    prob=sum1/Ns;

end