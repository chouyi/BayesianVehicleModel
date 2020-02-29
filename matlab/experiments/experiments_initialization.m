%% Initialization
dim=4;

% UAV model
if model == 1
    omega_max=0.30; 
    omega_min=-omega_max;
    para.max=omega_max;
    para.min=omega_min;

    [M,delimiterOut]=importdata('../../data/UAV_data.txt');
    [theta,delimiterOut]=importdata('../../data/theta_list.txt');

    if pred_step==5
        pre_com = csvread('../../outputs/ct-model-5-1.csv');
    elseif pred_step == 10
        pre_com = csvread('../../outputs/ct-model-10-1.csv');    
    elseif pred_step == 15
        pre_com = csvread('../../outputs/ct-model-15-1.csv');        
    end
    
    T=0.4;%discrete time interval dt
    para.num_p = 40;%number of omega particals
    np=1;% number of parameter;

    dOmega=[-0.045,0.045];
    
    dynamic.fun=@CTmodelDynamic; % dynamic functions
    dynamic.h=T;
    para.epsilon=0;%0.01;

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
    [M,delimiterOut]=importdata('../../data/car_data.txt');
    [theta,delimiterOut]=importdata('../../data/theta_list_Dubins.txt');
    pre_com = csvread('../../outputs/dubins-4-2.csv');
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

    dynamic.fun=@dubinsDynamic; % dynamic functions
    dynamic.h=dt;
    
    %---------------------------------------------
    idx.xx=1;   idx.xy=2;   idx.v=3;    idx.psi=4;

end
