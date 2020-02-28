clear all;
close all;

%% Input /
test_set_start=1;
test_set_end = 2;
Ns=20;% the number of testing samples of trajectory 

%% Load data
pred_step = 4;
pre_com = csvread('../../outputs/dubins-5-2.csv');
[M,delimiterOut]=importdata('../../data/car_data.txt');
[theta,delimiterOut]=importdata('../../data/theta_list_Dubins.txt');

%% params
dt=0.5;%discrete time interval dt
numInt = 10;%number of intervals/grind points on each parameter dimension
np=2;%paramer dim

maxtheta=[0.5 0.5];%[1 2];
mintheta=[-0.5 -0.5];%[-1 -2];

para.np=np;
para.num_p=numInt^np;
para.max= maxtheta;
para.min=mintheta; 
para.dpmax=[0.05 0.01];%[0.25 0.02];
para.dpmin=[-0.05 -0.01];%[-0.25 -0.02];
para.epsilon=0;%0.01;
para.theta=theta;
dim=4;

dynamic.fun=@dubinsDynamic; % dynamic functions
dynamic.h=dt;

NT=test_set_end- test_set_start;  

M_test = M(test_set_start:test_set_end+1+pred_step, :);
idx.xx=1;   idx.xy=2;   idx.v=3;    idx.psi=4;

para.prior_initial=(1/numInt)^np.*ones(1,numInt^np);%@(x) unifpdf(x,omega_min,omega_max);
para.prior=para.prior_initial;
para.epsilon=0.1;
%% Real-time part
result_report = zeros(1,3);
t= NT;

% Obstacle setting (at the end of a trajectory)
obs_center = [M_test(t+1+pred_step, idx.xx) M_test(t+1+pred_step, idx.xy)];
obs_zone = make_obstacle(obs_center, 0, 0);
result_report(1,1) = inpolygon(M_test(t+1+pred_step, idx.xx), M_test(t+1+pred_step, idx.xy), obs_zone.Vertices(:,1),obs_zone.Vertices(:,2));

tic;
%% Calc cov matrix
dynamic.cov=(1/100)*[[3, 0.1, 1, 1];[0.1, 3, 2, 1];[1, 2, 3, 0.1];[1, 1, 0.1, 3]];

x = zeros(dim, NT+1);
for t = 1:NT
    x(:,t) = [M_test(t,idx.xx); M_test(t,idx.xy); M_test(t,idx.v);M_test(t,idx.psi)];
    x(:,t+1) = [M_test(t+1,idx.xx); M_test(t+1,idx.xy); M_test(t+1,idx.v);M_test(t+1,idx.psi)];

    post=computePosterior_Dubinsmodel(para,x(:,t),x(:,t+1),dynamic);
    para.prior=(1-para.epsilon)*post+para.epsilon*para.prior_initial;

    % Only consider the last position
    if t== NT  
        sum1=0;
        for w_idx=1:para.num_p 
            w=theta(w_idx,:);
            pred_X=predict_ksteps_sampling_Dubins(Ns,w,para,x(:,t+1),dynamic,pred_step);
            temp_pc=Prob_collision(pred_X(1,:,pred_step),pred_X(2,:,pred_step),obs_zone);
            
            sum1=sum1 + temp_pc*post(w_idx);
            
            for i=1:size(pred_X,2)
                x1(1,:)=pred_X(1,i,:);
                y1(1,:)=pred_X(2,i,:);
                plot(x1,y1,'-b');hold on;
            end
        end
        result_report(1,2) = sum1;
    end
end
result_report(1,3) = toc;

patch(obs_zone.Vertices(:,1)', obs_zone.Vertices(:,2)', [0 0.4470 0.7410]);
hold on;
plot(M_test(2:NT+1+pred_step,idx.xx), M_test(2:NT+1+pred_step,idx.xy), '-*k');
hold on;    
plot(M_test(t+1+pred_step, idx.xx), M_test(t+1+pred_step, idx.xy), 'or');
hold on;
set(gca,'DataAspectRatio', [1 1 1]);
    