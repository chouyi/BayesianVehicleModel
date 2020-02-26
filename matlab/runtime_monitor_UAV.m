clear all
%%% Construct matrices for CT model
omega_max=0.18;
omega_min=-0.18;
T=0.4;%discrete time interval dt
para.num_p = 20;%number of omega particals
np=1;% number of parameter;

delta = (omega_max - omega_min)/(2*para.num_p);
CT.A=zeros(4,4,para.num_p);
for i=1:para.num_p
    W=omega_min + delta*(2*i-1); 
    WT=W*T;    
    CT.A(:,:,i)=[1 sin(WT)/W 0 -(1-cos(WT))/W; 0 cos(WT) 0 -sin(WT);...
                0 (1-cos(WT))/W 1 sin(WT)/W; 0 sin(WT) 0 cos(WT)];
end
% for i=1:para.num_p
%     W(i)=omega_min + delta*(2*i-1); 
% end
%%%-------------------
para.prior_initial=(1/para.num_p)*ones(1,para.num_p);%@(x) unifpdf(x,omega_min,omega_max);
para.prior=para.prior_initial;
para.epsilon=0;%0.01;
dim=4;
%p(omega(t)|x(t+1))=1/N*P(x(t+1))|omega)*Prior(omega(t))
%P(x(t+1))|omega)=N(Ax(t),sigma)


%%%%%%%%%%%%%
pre_com = csvread('../data/reach-sets-5-2.csv');
[M,delimiterOut]=importdata('../data/UAV_data.txt');
[theta,delimiterOut]=importdata('../data/theta_list.txt');
obs = importdata('../data/obstacles.txt');%%xmin,xmax,ymin,ymax
obs=obs(2,:);
CT.cov=zeros(dim,dim);

test_set_start=4280;%4200;
test_set_end = 4290;%4300;

training_set_start=test_set_start-1000;
training_set_end = test_set_start;

NT=test_set_end- test_set_start;  
NT_train=training_set_end- training_set_start;  

pre_kstep = 20;

% File name of data---------------------------
    M_test = M(test_set_start:test_set_end, :);
    M_training = M(training_set_start:training_set_end, :);
%---------------------------------------------

idx.xy = 3;         idx.xx = 4;         idx.xz = 5;
idx.vy = 6;         idx.vx = 7;         idx.vz = 8;
%% Calc cov matrix
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





%%%% obstacles %%%%

obs_coor_x = [obs(1) obs(1) obs(2) obs(2)];
obs_coor_y = [obs(3) obs(4) obs(4) obs(3)];
obs_zone = polyshape(obs_coor_x,obs_coor_y);
plot(obs_zone);
hold on;



x = zeros(dim, NT+1);

    for t = 1:NT
        
        x(:,t) = [M_test(t,idx.xx); M_test(t,idx.vx); M_test(t,idx.xy);M_test(t,idx.vy)];
        x(:,t+1) = [M_test(t+1,idx.xx); M_test(t+1,idx.vx); M_test(t+1,idx.xy);M_test(t+1,idx.vy)];
    
    %%compute posterior p(omega(t)|x(t+1))
        
        post=computePosterior_CTmodel(para,x(:,t),x(:,t+1),CT); 

        para.prior=post;%(1-para.epsilon)*post+para.epsilon*para.prior_initial;
        hold on
        plot(x(1,t+1),x(3,t+1), '*k'); 
        
        post_list(t,:)=post;
     
         

        %%predict prob of collision after k steps at kT
        sum1=0;
        for w_idx=1:para.num_p 
            w=theta(w_idx);
           
            collision = oracle(obs_zone,w, x(2,t+1), x(4,t+1), x(1,t+1), x(3,t+1), pre_com);
            if collision==1
                sum1=sum1 + post(w_idx);
            end
            prob_c(t)=sum1;
        end
%         if (mod(t,pre_kstep)==0)    
%         pred_X=predict_ksteps(para,x(:,t+1),CT,pre_kstep);
%         %para.prior=para.prior_initial;
%         hold on;
%         for i=1:size(pred_X,2)
% 
%             x1(1,:)=pred_X(1,i,:);
%             y1(1,:)=pred_X(3,i,:);
% 
%            plot(x1,y1,'-b')
%             %plot(pred_X(1,i,:),pred_X(3,i,:),'-ob')
%            
%         end
%        end%if

    end
    

    
    