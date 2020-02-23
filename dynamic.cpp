#include "dynamic.h"
vector<double> CTmodel_dynamic(vector<double> cur_x,vector<double> theta,double dt) {
    //h: time step size
    int sdim=cur_x.size();
    
    vector<double> next_x(sdim);
    double W=theta[0];
    double WT=W*dt;
    if(W==0){
        cout<<"error:omega=0"<<endl;
    }
    next_x[0]=cur_x[0] + sin(WT)/W*cur_x[1] + (1-cos(WT))/W*cur_x[3];
    next_x[1]= cos(WT)*cur_x[1] -sin(WT)*cur_x[3];
    next_x[2]=(1-cos(WT))/W*cur_x[1] + cur_x[2] + sin(WT)/W*cur_x[3];
    next_x[3]= sin(WT)*cur_x[1] + cos(WT)*cur_x[3];
    //cout<<next_x[0]<<endl;
    //    WT=W*T;
    //    CT.A(:,:,i)=[1 sin(WT)/W 0 (1-cos(WT))/W; 0 cos(WT) 0 -sin(WT);...
    //                 0 (1-cos(WT))/W 1 sin(WT)/W; 0 sin(WT) 0 cos(WT)];
    
    return next_x;
    
}
