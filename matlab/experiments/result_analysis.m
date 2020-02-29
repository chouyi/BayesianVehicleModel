clear all;
threshold = 0.950;
%UAV
pred_step_set = [5];
update_step_set = [20];
%Dubins car
% pred_step_set_dubins = [4];
% update_step_set_dubins = [0 10 20 50];

outputname_set = [string('./report/UAV_reachable_set-') string('./report/UAV_sampling-')...
              string('./report/dubins_reachable_set-') string('./report/dubins_sampling-')];

total_report = zeros(15,8);
total_report_UAV_reachable_set = zeros(15,2);
total_report_UAV_sampling = zeros(15,2);
total_report_dubins_reachable_set = zeros(15,2);
total_report_dubins_sampling = zeros(15,2);
% 
% for k = 1:1
%     outputname = outputname_set(k);
% for j = 1:size(pred_step_set,2)
%     for l = 1:size(update_step_set,2)
%         update_step = update_step_set(l);  
%         pred_step = pred_step_set(j); 
%         filename = strcat(outputname,string(update_step), '-', string(pred_step),'result_report.csv');
%         
%         result = csvread(filename);
%         computation_time = mean(result(:, 5));
%         correct = 0;
%         wrong = 0;
% 
%         for i = 1:size(result,1)
%             if result(i,2) == 1 && result(i,3) >= threshold 
%                 correct = correct + 1;
%             elseif result(i,2) == 1 && result(i,3) < threshold 
%                 wrong = wrong + 1;
%             end
%         end
% 
%         rate_correct = correct / (correct+wrong);
%         rate_wrong = wrong / (correct+wrong);
%         
%         total_report(l+(j-1)*size(update_step_set,2), k+2*(k-1)) = rate_correct;
%         total_report(l+(j-1)*size(update_step_set,2), k+2*(k-1)+1) = rate_wrong;
%     end
% end
%     total_report_UAV_reachable_set = total_report(:,1:2);
%     total_report_UAV_sampling = total_report(:,3:4);
% end

% 
% for k = 3:4
%     outputname = outputname_set(k);
% for j = 1:size(pred_step_set_dubins,2)
%     for l = 1:size(update_step_set_dubins,2)
%         update_step = update_step_set_dubins(l);  
%         pred_step = pred_step_set_dubins(j); 
%         filename = strcat(outputname,string(update_step), '-', string(pred_step),'result_report.csv');
%         
%         result = csvread(filename);
%         computation_time = mean(result(:, 5));
%         correct = 0;
%         wrong = 0;
% 
%         for i = 1:size(result,1)
%             if result(i,2) == 1 && result(i,3) >= threshold 
%                 correct = correct + 1;
%             elseif result(i,2) == 1 && result(i,3) < threshold 
%                 wrong = wrong + 1;
%             end
%         end
% 
%         rate_correct = correct / (correct+wrong);
%         rate_wrong = wrong / (correct+wrong);
%         
%         total_report(l+(j-1)*size(update_step_set,2), k+2*(k-1)) = rate_correct;
%         total_report(l+(j-1)*size(update_step_set,2), k+2*(k-1)+1) = rate_wrong;
%     end
% end
%     total_report_dubins_reachable_set = total_report(:,5:6);
%     total_report_dubins_sampling = total_report(:,7:8);
% end
% 
% 
for k = 1:1
    outputname = outputname_set(k);
for j = 1:size(pred_step_set,2)
    for l = 1:size(update_step_set,2)
        update_step = update_step_set(l);  
        pred_step = pred_step_set(j); 
        filename = strcat(outputname,string(update_step), '-', string(pred_step),'result_report.csv');
        
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
        
        total_report(l+(j-1)*size(update_step_set,2), k+2*(k-1)) = rate_correct;
        total_report(l+(j-1)*size(update_step_set,2), k+2*(k-1)+1) = rate_wrong;
    end
end
    total_report_UAV_reachable_set = total_report(:,1:2);
    total_report_UAV_sampling = total_report(:,3:4);
end