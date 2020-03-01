clear all;
close all;

%% Input 
test_set_start=2721;%4280;%4200;
test_set_end = 2731;%4290;%4300;

%% Load data
pred_step = 10;
pre_com = csvread('../../outputs/ct-model-5-2.csv');
[M,delimiterOut]=importdata('../../data/UAV_data.txt');
[theta,delimiterOut]=importdata('../../data/theta_list.txt');

%% params
omega_max=0.3; 
omega_min=-0.3;

T=0.4;%discrete time interval dt
para.num_p = 40;%number of omega particals
para.max=omega_max;
para.min=omega_min;
np=1;% number of parameter;
Ns=100;% the number of testing samples of trajectory  
dOmega=[-0.045,0.045];
delta = (omega_max - omega_min)/(2*para.num_p);
dynamic.fun=@CTmodelDynamic; % dynamic functions
dynamic.h=T;
para.epsilon=0.01;
dim=4;

%% Construct CT model
delta = (omega_max - omega_min)/(2*para.num_p);
CT.A=zeros(4,4,para.num_p);
for i=1:para.num_p
    W=omega_min + delta*(2*i-1); 
    WT=W*T;    
    CT.A(:,:,i)=[1 sin(WT)/W 0 -(1-cos(WT))/W; 0 cos(WT) 0 -sin(WT);...
                0 (1-cos(WT))/W 1 sin(WT)/W; 0 sin(WT) 0 cos(WT)];
end

training_set_start=test_set_start-1000;
training_set_end = test_set_start;

NT=test_set_end- test_set_start;  
NT_train=training_set_end- training_set_start;  

M_test = M(test_set_start:test_set_end+1+pred_step, :);
M_training = M(training_set_start:training_set_end+pred_step, :);
idx.xy = 3;         idx.xx = 4;         idx.xz = 5;
idx.vy = 6;         idx.vx = 7;         idx.vz = 8;

para.prior_initial=(1/para.num_p)*ones(1,para.num_p);%@(x) unifpdf(x,omega_min,omega_max);
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
CT.cov=zeros(dim,dim);
for i = 1:NT_train-1
    W_train(i) = atan2(M_training(i+1, idx.vy),M_training(i+1, idx.vx))...
            - atan2(M_training(i, idx.vy)  ,M_training(i, idx.vx));
    if W_train(i) > pi
        W_train(i) = W_train(i) - ceil(W_train(i)/pi)*pi;
    elseif W_train(i) < -pi
        W_train(i) = W_train(i) + ceil(W_train(i)/pi)*pi;
    end
   
    WT=W_train(i)*T;
    
    CT_model = [1 sin(WT)/W_train(i) 0 -(1-cos(WT))/W_train(i); 
                0 cos(WT) 0 -sin(WT);...
                0 (1-cos(WT))/W_train(i) 1 sin(WT)/W_train(i); 
                0 sin(WT) 0 cos(WT)];
            
    error_train(:,i) = [M_training(i+1, idx.xx);M_training(i+1, idx.vx);...
                        M_training(i+1, idx.xy);M_training(i+1, idx.vy)]...
                        -CT_model * [M_training(i, idx.xx);...
                                    M_training(i, idx.vx);...
                                    M_training(i, idx.xy);...
                                    M_training(i, idx.vy)]; 
end
CT.cov = cov(error_train');

x = zeros(dim, NT+1);
for t = 1:NT
    x(:,t) = [M_test(t,idx.xx); M_test(t,idx.vx); M_test(t,idx.xy);M_test(t,idx.vy)];
    x(:,t+1) = [M_test(t+1,idx.xx); M_test(t+1,idx.vx); M_test(t+1,idx.xy);M_test(t+1,idx.vy)];

    post=computePosterior_CTmodel(para,x(:,t),x(:,t+1),CT); 
    para.prior=(1-para.epsilon)*post+para.epsilon*para.prior_initial;
    
    % Only consider the last position
    if t== NT  
        sum1=0;
        post_theta = [post;theta']';
        post_theta = sortrows(post_theta,1);
        
        for w_idx=31:40%para.num_p 
            w=post_theta(w_idx,2);
            
            pred_X=predict_ksteps_sampling_UAV(Ns,dOmega,w,para,x(:,t+1),dynamic,pred_step);
            temp_pc=Prob_collision(pred_X(1,:,pred_step),pred_X(3,:,pred_step),obs_zone);
            
            sum1=sum1 + temp_pc*post_theta(w_idx,1);
            for i=1:size(pred_X,2)
                x1(1,:)=pred_X(1,i,:);
                y1(1,:)=pred_X(3,i,:);
                plot(x1,y1,'-b');
                hold on;
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
    