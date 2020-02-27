%% Initialization
dim=4;
% mkdir results_figures/reachable_set/collision 0-5
% mkdir results_figures/reachable_set/collision 0-10
% mkdir results_figures/reachable_set/collision 0-15
% mkdir results_figures/reachable_set/collision 10-5
% mkdir results_figures/reachable_set/collision 10-10
% mkdir results_figures/reachable_set/collision 10-15
% mkdir results_figures/reachable_set/collision 20-5
% mkdir results_figures/reachable_set/collision 20-10
% mkdir results_figures/reachable_set/collision 20-15
% mkdir results_figures/reachable_set/collision 50-5
% mkdir results_figures/reachable_set/collision 50-10
% mkdir results_figures/reachable_set/collision 50-15
% mkdir results_figures/reachable_set/collision 100-5
% mkdir results_figures/reachable_set/collision 100-10
% mkdir results_figures/reachable_set/collision 100-15

mkdir results_figures/reachable_set/no_collision 0-5
mkdir results_figures/reachable_set/no_collision 0-10
mkdir results_figures/reachable_set/no_collision 0-15
mkdir results_figures/reachable_set/no_collision 10-5
mkdir results_figures/reachable_set/no_collision 10-10
mkdir results_figures/reachable_set/no_collision 10-15
mkdir results_figures/reachable_set/no_collision 20-5
mkdir results_figures/reachable_set/no_collision 20-10
mkdir results_figures/reachable_set/no_collision 20-15
mkdir results_figures/reachable_set/no_collision 50-5
mkdir results_figures/reachable_set/no_collision 50-10
mkdir results_figures/reachable_set/no_collision 50-15
mkdir results_figures/reachable_set/no_collision 100-5
mkdir results_figures/reachable_set/no_collision 100-10
mkdir results_figures/reachable_set/no_collision 100-15
% 
% mkdir results_figures/sampling/collision 0-5
% mkdir results_figures/sampling/collision 0-10
% mkdir results_figures/sampling/collision 0-15
% mkdir results_figures/sampling/collision 10-5
% mkdir results_figures/sampling/collision 10-10
% mkdir results_figures/sampling/collision 10-15
% mkdir results_figures/sampling/collision 20-5
% mkdir results_figures/sampling/collision 20-10
% mkdir results_figures/sampling/collision 20-15
% mkdir results_figures/sampling/collision 50-5
% mkdir results_figures/sampling/collision 50-10
% mkdir results_figures/sampling/collision 50-15
% mkdir results_figures/sampling/collision 100-5
% mkdir results_figures/sampling/collision 100-10
% mkdir results_figures/sampling/collision 100-15
% 
% mkdir results_figures/sampling/no_collision 0-5
% mkdir results_figures/sampling/no_collision 0-10
% mkdir results_figures/sampling/no_collision 0-15
% mkdir results_figures/sampling/no_collision 10-5
% mkdir results_figures/sampling/no_collision 10-10
% mkdir results_figures/sampling/no_collision 10-15
% mkdir results_figures/sampling/no_collision 20-5
% mkdir results_figures/sampling/no_collision 20-10
% mkdir results_figures/sampling/no_collision 20-15
% mkdir results_figures/sampling/no_collision 50-5
% mkdir results_figures/sampling/no_collision 50-10
% mkdir results_figures/sampling/no_collision 50-15
% mkdir results_figures/sampling/no_collision 100-5
% mkdir results_figures/sampling/no_collision 100-10
% mkdir results_figures/sampling/no_collision 100-15
% UAV model
if model == 1
    omega_max=0.18; 
    omega_min=-omega_max;
    para.max=omega_max;
    para.min=omega_min;

    [M,delimiterOut]=importdata('../data/UAV_data.txt');
    [theta,delimiterOut]=importdata('../data/theta_list.txt');

    if pred_step==5
        pre_com = csvread('../outputs/ct-model-5-2.csv');
    elseif pred_step == 10
        pre_com = csvread('../outputs/ct-model-10-2.csv');    
    elseif pred_step == 15
        pre_com = csvread('../outputs/ct-model-15-2.csv');        
    end
    
    T=0.4;%discrete time interval dt
    para.num_p = 20;%number of omega particals
    
    %---------------------------------------------
    idx.xy = 3;         idx.xx = 4;         idx.xz = 5;
    idx.vy = 6;         idx.vx = 7;         idx.vz = 8;

    %% Construct CT model
    delta = (omega_max - omega_min)/(2*para.num_p);
    CT.A=zeros(4,4,para.num_p);
    for i=1:para.num_p
        W=omega_min + delta*(2*i-1); 
        WT=W*T;    
        CT.A(:,:,i)=[1 sin(WT)/W 0 -(1-cos(WT))/W; 0 cos(WT) 0 -sin(WT);...
                    0 (1-cos(WT))/W 1 sin(WT)/W; 0 sin(WT) 0 cos(WT)];
    end
    
% Dubins car model
elseif model ==2
    [M,delimiterOut]=importdata('../data/car_data.txt');
    [theta,delimiterOut]=importdata('../data/theta_list_Dubins.txt');
    pre_com = csvread('../outputs/output.csv');
    
    T=0.5;%discrete time interval dt
    numInt = 10;%number of intervals/grind points on each parameter dimension
    np=2;%paramer dim
    para.num_p=numInt^np;
    para.theta=theta;
    
    test_set_start=10;
    test_set_end = 11;

    dynamic.fun=@dubinsDynamic; % dynamic functions
    dynamic.h=T;

    % File name of data---------------------------
    M_test = M(test_set_start:test_set_end+pred_step, :);
    %---------------------------------------------
    idx.xx=1;   idx.xy=2;   idx.v=3;    idx.psi=4;

    para.prior_initial=(1/numInt)^np.*ones(1,numInt^np);%@(x) unifpdf(x,omega_min,omega_max);
    para.prior=para.prior_initial;
    NT=test_set_end- test_set_start;  
end
