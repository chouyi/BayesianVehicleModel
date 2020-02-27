clear all
dt=0.5;%discrete time interval dt
numInt = 10;%number of intervals/grind points on each parameter dimension
np=2;%paramer dim
Ns=20;% the number of testing samples of trajectory 
para.np=np;
para.num_p=numInt^np;
pre_kstep = 4;
test_set_start=1;
test_set_end = 20;
maxtheta=[0.5 0.5];%[1 2];
mintheta=[-0.5 -0.5];%[-1 -2];
para.max= maxtheta;
para.min=mintheta; 
para.dpmax=[0.05 0.01];%[0.25 0.02];
para.dpmin=[-0.05 -0.01];%[-0.25 -0.02];

dynamic.fun=@dubinsDynamic; % dynamic functions
dynamic.h=dt;
% training_set_start=1;
% training_set_end = 50;

NT=test_set_end- test_set_start;  
% NT_train=training_set_end- training_set_start;  

filename = '../data/car_data.txt';
[M,delimiterOut]=importdata(filename);
% filename = 'logPosterior_list.txt';
% [logP,delimiterOut]=importdata(filename);
filename = '../data/theta_list_Dubins.txt';
[theta,delimiterOut]=importdata(filename);
pre_com = csvread('../outputs/dubins-5-2.csv');

obs = [-25 -5 15 25];

para.theta=theta;

idx.xx=1;
idx.xy=2;
idx.v=3;
idx.psi=4;

para.prior_initial=(1/numInt)^np.*ones(1,numInt^np);%@(x) unifpdf(x,omega_min,omega_max);

para.prior=para.prior_initial;
%para.epsilon=0.005;
dim=4;

dynamic.cov=(1/100)*[[3, 0.1, 1, 1];[0.1, 3, 2, 1];[1, 2, 3, 0.1];[1, 1, 0.1, 3]];
%%%% obstacles %%%%

obs_coor_x = [obs(1) obs(1) obs(2) obs(2)];
obs_coor_y = [obs(3) obs(4) obs(4) obs(3)];
obs_zone = polyshape(obs_coor_x,obs_coor_y);
plot(obs_zone);
hold on;

M_test = M(test_set_start:test_set_end, :);
% M_training = M(training_set_start:training_set_end, :);

x = zeros(dim, NT+1);

    for t = 1:NT
        x(:,t) = [M_test(t,idx.xx); M_test(t,idx.xy); M_test(t,idx.v);M_test(t,idx.psi)];
        x(:,t+1) = [M_test(t+1,idx.xx); M_test(t+1,idx.xy); M_test(t+1,idx.v);M_test(t+1,idx.psi)];
        
        post=computePosterior_Dubinsmodel(para,x(:,t),x(:,t+1),dynamic);


        para.prior=post;%(1-para.epsilon)*post+para.epsilon*para.prior_initial;
        hold on
        plot(x(1,t+1),x(2,t+1), '*k'); 
       % plot(para.prior)
      %%predict prob of collision after k steps at kT
%         sum1=0;
%         for w_idx=1:para.num_p  
%             w=theta(w_idx,:);
% 
%             collision = oracle_Dubins(obs_zone,w, x(3,t+1), x(4,t+1), x(1,t+1), x(2,t+1), pre_com);
%             if collision==1
%                 sum1=sum1 + post(w_idx);
%             end
%             prob_c(t)=sum1;
%         end
        %%predict k steps at kT
        sum1=0;
        for w_idx=1:para.num_p 
            w=theta(w_idx,:);

            pred_X=predict_ksteps_sampling_Dubins(Ns,w,para,x(:,t+1),dynamic,pre_kstep);
            temp_pc=Prob_collision(pred_X(1,:,pre_kstep),pred_X(2,:,pre_kstep),obs);
            
            sum1=sum1 + temp_pc*post(w_idx);
            
            hold on;
            for i=1:size(pred_X,2)

                x1(1,:)=pred_X(1,i,:);
                y1(1,:)=pred_X(2,i,:);

                plot(x1,y1,'-b')
            %plot(pred_X(1,i,:),pred_X(3,i,:),'-ob')
           
            end
        end%w
        prob_c(t)=sum1;%%the collision probablilty over time
    end
