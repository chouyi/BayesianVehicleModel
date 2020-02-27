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

%% Setting
pred_step_set = [4];
update_step_set = [0 10 20 50];
ptime = 2; % pause time to see figures before closing it

model = 2; %dubins=2
num_test = 2;

obs_size = 4;
obs_dist_min = 0;
obs_dist_max = 0;

output_filename = strcat('./result_report.csv');

%% 
for j = 1:size(pred_step_set,2)
    for l = 1:size(update_step_set,2)
    update_step = update_step_set(l);  
    pred_step = pred_step_set(j);  
    
    NT = update_step+1;
    experiments_initialization;

    %% Real-time part
    % Simulation
    result_report = zeros(1,3);
    for k = 1:num_test
        para.prior_initial=(1/numInt)^np.*ones(1,numInt^np);
        para.prior=para.prior_initial;

        % test point
        test_set_start = randi([0 400],1,1);
        result_report(k,1) = test_set_start;
        M_test = M(test_set_start:test_set_start+1+pred_step+NT, :);

        t= NT; % to make an obstacle at the end of trajectory
        % Obstacle setting
        obs_center = [M_test(t+1+pred_step, idx.xx) M_test(t+1+pred_step, idx.xy)];
        obs_zone = make_obstacle(obs_center, obs_size, obs_dist_min, obs_dist_max);
        result_report(k,2) = inpolygon(M_test(t+1+pred_step, idx.xx), M_test(t+1+pred_step, idx.xy), obs_zone.Vertices(:,1),obs_zone.Vertices(:,2));

        % Record a distance from obstacle
        result_report(k,4) = min(sqrt((obs_zone.Vertices(:,1)-M_test(t+1+pred_step, idx.xx)).^2+...
                            (obs_zone.Vertices(:,2)-M_test(t+1+pred_step, idx.xy)).^2));

        tic;
        dynamic.cov=(1/100)*[[3, 0.1, 1, 1];[0.1, 3, 2, 1];[1, 2, 3, 0.1];[1, 1, 0.1, 3]];
        x = zeros(dim, NT+1);

        for t = 1:NT
            x(:,t) = [M_test(t,idx.xx); M_test(t,idx.xy); M_test(t,idx.v);M_test(t,idx.psi)];
            x(:,t+1) = [M_test(t+1,idx.xx); M_test(t+1,idx.xy); M_test(t+1,idx.v);M_test(t+1,idx.psi)];

            post=computePosterior_Dubinsmodel(para,x(:,t),x(:,t+1),dynamic);
            para.prior=post;%(1-para.epsilon)*post+para.epsilon*para.prior_initial;

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
                result_report(k,3)  = sum1;
            end
        end

        result_report(k,5) = toc;

        patch(obs_zone.Vertices(:,1)', obs_zone.Vertices(:,2)', [0 0.4470 0.7410]);
        hold on;
        plot(M_test(2:NT+1+pred_step,idx.xx), M_test(2:NT+1+pred_step,idx.xy), '-*k');
        hold on;    
        plot(M_test(t+1+pred_step, idx.xx), M_test(t+1+pred_step, idx.xy), 'or');
        hold on;
    %     filename = strcat('./figures/',string(update_step), '-', string(pred_step), '-', string(k), '.png');
    %     filename_fig = strcat('./figures/',string(update_step), '-', string(pred_step), '-', string(k), '.fig');
    %     saveas(gcf,filename)
    %     saveas(gcf,filename_fig)
        set(gca,'DataAspectRatio', [1 1 1]);
        pause(ptime);
        close all;
    end
        writematrix(result_report,output_filename)
    
    end
end
