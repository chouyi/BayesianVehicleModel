clear all;
threshold = 0.950;

total_report_collision = zeros(15,2);
total_report_no_collision = zeros(15,2);
% 1. Different update time to allow some time to update posterior
% 2. Different steps - collision
pred_step_set = [5 10 15];
update_step_set = [0 10 20 50 100];

for j = 1:size(pred_step_set,2)
    for l = 1:size(update_step_set,2)
        update_step = update_step_set(l);  
        pred_step = pred_step_set(j); 
        filename = strcat('./results_figures/reachable_set/collision/',string(update_step), '-', string(pred_step), '/result_report2.csv');
        result = csvread(filename);
        computation_time = mean(result(:, 5));
        correct = 0;
        wrong = 0;

        for i = 1:size(result,1)
            if result(i,2) == 1 && result(i,3) >= threshold 
                correct = correct + 1;
            elseif result(i,2) == 1 && result(i,3) < threshold 
                wrong = wrong + 1;
            end
        end

        rate_correct = correct / (correct+wrong);
        rate_wrong = wrong / (correct+wrong);
        
        total_report_collision(l+(j-1)*size(update_step_set,2), 1) = rate_correct;
        total_report_collision(l+(j-1)*size(update_step_set,2), 2) = rate_wrong;
    end
end
% 3. Different steps ? no collision (close enough, random distance between 1-10)
clear result rate_correct rate_wrong
for j = 1:size(pred_step_set,2)
    for l = 1:size(update_step_set,2)
        update_step = update_step_set(l);  
        pred_step = pred_step_set(j); 
        filename = strcat('./results_figures/reachable_set/no_collision/',string(update_step), '-', string(pred_step), '/result_report3.csv');
        result = csvread(filename);
        computation_time = mean(result(:, 5));
        correct = 0;
        wrong = 0;

        for i = 1:size(result,1)
            if result(i,2) == 0 && result(i,3) < threshold 
                correct = correct + 1;
            elseif result(i,2) == 0 && result(i,3) >= threshold 
                wrong = wrong + 1;
            end
        end

        rate_correct = correct / (correct+wrong);
        rate_wrong = wrong / (correct+wrong);
        
        total_report_no_collision(l+(j-1)*size(update_step_set,2), 1) = rate_correct;
        total_report_no_collision(l+(j-1)*size(update_step_set,2), 2) = rate_wrong;
    end
end