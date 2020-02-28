clear all;
close all;

%% Input 
update_step = 10;
test_set_start=1;

%% Load data
pred_step = 4;
pre_com = csvread('../../outputs/dubins-4-2.csv');
[M,delimiterOut]=importdata('../../data/car_data.txt');
[theta,delimiterOut]=importdata('../../data/theta_list_Dubins.txt');

%% params
T=0.5;%discrete time interval dt
numInt = 10;%number of intervals/grind points on each parameter dimension
np=2;%paramer dim
para.num_p=numInt^np;
dim=4;
para.theta=theta;
dynamic.fun=@dubinsDynamic; % dynamic functions
dynamic.h=T;
NT = update_step+1;

M_test = M(test_set_start:test_set_start+1+pred_step+NT, :);
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
        post_theta = [post;theta']';
        post_theta = sortrows(post_theta,1);
        for w_idx=1:para.num_p  
            w=post_theta(w_idx,2:end);
            [collision, reach_set] = oracle_Dubins(obs_zone, w(1), w(2), x(3,t+1), x(4,t+1), x(1,t+1), x(2,t+1), pre_com);
            rgb = [1-post_theta(w_idx,1) 1-post_theta(w_idx,1) 1-post_theta(w_idx,1)];

            if collision==1
                sum1=sum1 + post_theta(w_idx,1);
                rgb = [1-post_theta(w_idx,1) 0 0];
            end

            if post_theta(w_idx,1) < 0.01
                continue;
            end
            patch(reach_set.Vertices(:,1)', reach_set.Vertices(:,2)', rgb, 'EdgeColor','none');
            hold on;
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
    