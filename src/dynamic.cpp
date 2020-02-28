#include "dynamic.h"
vector<double> CTmodel_dynamic(vector<double> cur_x,vector<double> theta,double dt) {
    //dt: time step size
    int sdim=cur_x.size();
    
    vector<double> next_x(sdim);
    double W=theta[0];
    double WT=W*dt;
    if(W==0){
        next_x[0]= cur_x[0] + cur_x[1] ;
        next_x[1]= cos(WT)*cur_x[1] -sin(WT)*cur_x[3];
        next_x[2]= cur_x[2] + cur_x[3];
        next_x[3]= sin(WT)*cur_x[1] + cos(WT)*cur_x[3];
    }else{
        next_x[0]=cur_x[0] + sin(WT)/W*cur_x[1] - (1-cos(WT))/W*cur_x[3];
        next_x[1]= cos(WT)*cur_x[1] -sin(WT)*cur_x[3];
        next_x[2]=(1-cos(WT))/W*cur_x[1] + cur_x[2] + sin(WT)/W*cur_x[3];
        next_x[3]= sin(WT)*cur_x[1] + cos(WT)*cur_x[3];
    }
    //cout<<next_x[0]<<endl;
    //    WT=W*T;
    //    CT.A(:,:,i)=[1 sin(WT)/W 0 -(1-cos(WT))/W; 0 cos(WT) 0 -sin(WT);...
    //                 0 (1-cos(WT))/W 1 sin(WT)/W; 0 sin(WT) 0 cos(WT)];
    
    return next_x;
    
}

vector<double> DubinsModel_dynamic(vector<double> x,vector<double> theta,double dt) {
    //h: time step size
    double h=dt;
    int sdim= x.size();
    
    vector<double> next_x(sdim);
    double u_v=theta[0];
    double u_psi=theta[1];

    next_x[0]=x[0]+h*x[2]*cos(x[3])+(h*h/2)*u_v*cos(x[3]);
    next_x[1]=x[1]+h*x[2]*sin(x[3])+(h*h/2)*u_v*sin(x[3]);
    next_x[2]=x[2]+h*u_v;
    next_x[3]=x[3]+h*u_psi;

    return next_x;
    
}
