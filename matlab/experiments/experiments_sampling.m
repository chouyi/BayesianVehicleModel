clear;
close all;
%%
% Test setting
pred_step_set = [5 10 15];
update_step_set = [0 10 20 50 100];

for j = 1:size(pred_step_set,2)
    for l = 1:size(update_step_set,2)
update_step = update_step_set(l);  
pred_step = pred_step_set(j);  
model = 1; %UAV=1
num_test = 500;

% test_set_start=4270;%4200;
% test_set_end = 4271;%4300;

experiments_initialization;
Ns=100;% the number of testing samples of trajectory  
dOmega=[-0.045,0.045];
dynamic.fun=@CTmodelDynamic; % dynamic functions
dynamic.h=T;
para.epsilon=0;%0.01;

%% Real-time part

% Simulation

NT = update_step+1;

result_report = zeros(1,3);
for k = 1:num_test
    para.prior_initial=(1/para.num_p)*ones(1,para.num_p);%@(x) unifpdf(x,omega_min,omega_max);
    para.prior=para.prior_initial;
    
    test_set_start = randi([3200 5200],1,1);

    result_report(k,1) = test_set_start;
    M_test = M(test_set_start:test_set_start+1+pred_step+NT, :);
    
    training_set_start=test_set_start-1000;
    training_set_end = test_set_start;
    NT_train=training_set_end- training_set_start;  

    % Calc cov matrix
    CT.cov=zeros(dim,dim);
    M_training = M(training_set_start:training_set_end, :);

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
    
t= NT;
% Obstacle setting
obs_center = [M_test(t+1+pred_step, idx.xx) M_test(t+1+pred_step, idx.xy)];
obs_zone = make_obstacle(obs_center);
result_report(k,2) = inpolygon(M_test(t+1+pred_step, idx.xx), M_test(t+1+pred_step, idx.xy), obs_zone.Vertices(:,1),obs_zone.Vertices(:,2));

% Distance from obstacle
result_report(k,4) = min(sqrt((obs_zone.Vertices(:,1)-M_test(t+1+pred_step, idx.xx)).^2+...
                    (obs_zone.Vertices(:,2)-M_test(t+1+pred_step, idx.xy)).^2));
                        
tic;
x = zeros(dim, NT+1);
for t = 1:NT
    x(:,t) = [M_test(t,idx.xx); M_test(t,idx.vx); M_test(t,idx.xy);M_test(t,idx.vy)];
    x(:,t+1) = [M_test(t+1,idx.xx); M_test(t+1,idx.vx); M_test(t+1,idx.xy);M_test(t+1,idx.vy)];
    
    %%compute posterior p(omega(t)|x(t+1))
    post=computePosterior_CTmodel(para,x(:,t),x(:,t+1),CT); 
    para.prior=post;%(1-para.epsilon)*post+para.epsilon*para.prior_initial;

    if t== NT
        %%predict prob of collision after k steps at kT
        sum1=0;
        post_theta = [post;theta']';
        post_theta = sortrows(post_theta,1);
        for w_idx=1:para.num_p 
            w=post_theta(w_idx,2);
            pred_X=predict_ksteps_sampling_UAV(Ns,dOmega,w,para,x(:,t+1),dynamic,pred_step);
            temp_pc=Prob_collision(pred_X(:,:,pred_step),obs_zone);
            sum1 = sum1 + temp_pc*post_theta(w_idx,2);
            
            rgb = [1-sum1 1-sum1 1-sum1];

            if sum1>=0.95
                rgb = [1-post_theta(w_idx,1) 0 0];
            end

            if sum1 < 0.05
                continue;
            end
%             patch(reach_set.Vertices(:,1)', reach_set.Vertices(:,2)', rgb, 'EdgeColor','none');
%             hold on;
%             for i=1:size(pred_X,2)
% 
%                 x1(1,:)=pred_X(1,i,:);
%                 y1(1,:)=pred_X(3,i,:);

%                 plot(x1,y1,'-b');
%                 hold on;
%             end
        end
        result_report(k,3) = sum1;
    end
end

result_report(k,5) = toc;

% patch(obs_zone.Vertices(:,1)', obs_zone.Vertices(:,2)', [0 0.4470 0.7410]);
% hold on;
% plot(M_test(2:NT+1+pred_step,idx.xx), M_test(2:NT+1+pred_step,idx.xy), '-*k');
% hold on;    
% plot(M_test(t+1+pred_step, idx.xx), M_test(t+1+pred_step, idx.xy), 'or');
% hold on;
% filename = strcat('./results_figures/',string(update_step), '-', string(pred_step), '/', string(k), '.png');
% filename_fig = strcat('./results_figures/',string(update_step), '-', string(pred_step), '/', string(k), '.fig');
% saveas(gcf,filename)
% saveas(gcf,filename_fig)
% set(gca,'DataAspectRatio', [1 1 1]);
% close all;
end
filename = strcat('./results_figures/',string(update_step), '-', string(pred_step), '/result_report4.csv');
writematrix(result_report,filename)
    end
end
