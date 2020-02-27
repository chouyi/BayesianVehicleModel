clear;
close all;
%%
% Test: Collision detection -> accuracy including false alarm rate, time performance
% 1-1: 5 steps, collision, close but not collision
% 1-2: 10 steps, collision, close but not collision
% 1-3: 15 steps, collision, close but not collision
% 1-4: Combined test: use three steps together and show how they support
% each other by removing blind spot.

% Compare with sampling
% Different update time-compare result

% Close apporach test: various distances between the trajectory and obstacle
% Collision test: middle of an obstacle, edge of an obstacle

%%
% Test setting
pred_step_set = [5 10 15];
update_step_set = [0 10 20 50 100];

for j = 1:size(pred_step_set,2)
    for l = 1:size(update_step_set,2)
update_step = update_step_set(l);  
pred_step = pred_step_set(j);  
model = 1; %UAV=1
num_test = 100;

% test_set_start=4270;%4200;
% test_set_end = 4271;%4300;

experiments_initialization;

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
            [collision, reach_set] = oracle(obs_zone, w, x(2,t+1), x(4,t+1), x(1,t+1), x(3,t+1), pre_com);
            rgb = [1-post_theta(w_idx,1) 1-post_theta(w_idx,1) 1-post_theta(w_idx,1)];

            if collision==1
                sum1=sum1 + post_theta(w_idx,1);
                rgb = [1-post_theta(w_idx,1) 0 0];
            end

            if post_theta(w_idx,1) < 0.05
                continue;
            end
            patch(reach_set.Vertices(:,1)', reach_set.Vertices(:,2)', rgb, 'EdgeColor','none');
            hold on;
        end
        result_report(k,3) = sum1;
    end
end

result_report(k,5) = toc;

patch(obs_zone.Vertices(:,1)', obs_zone.Vertices(:,2)', [0 0.4470 0.7410]);
hold on;
plot(M_test(2:NT+1+pred_step,idx.xx), M_test(2:NT+1+pred_step,idx.xy), '-*k');
hold on;    
plot(M_test(t+1+pred_step, idx.xx), M_test(t+1+pred_step, idx.xy), 'or');
hold on;
filename = strcat('./results_figures/reachable_set/no_collision/',string(update_step), '-', string(pred_step), '/', string(k), '.png');
filename_fig = strcat('./results_figures/reachable_set/no_collision/',string(update_step), '-', string(pred_step), '/', string(k), '.fig');
saveas(gcf,filename)
saveas(gcf,filename_fig)
set(gca,'DataAspectRatio', [1 1 1]);
close all;
end
filename = strcat('./results_figures/reachable_set/no_collision/',string(update_step), '-', string(pred_step), '/result_report3.csv');
writematrix(result_report,filename)
    end
end
