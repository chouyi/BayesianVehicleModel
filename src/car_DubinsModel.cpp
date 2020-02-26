
#include "inference.h"

#include <sstream>      // std::istringstream
#include <fstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include <numeric>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <algorithm>


using namespace std;

//const string ctModelDataFile = "../data/car_data.txt";
int main(){
    vector<double> maxtheta,mintheta,delta_p;
    maxtheta.push_back(1);
    maxtheta.push_back(2);
    mintheta.push_back(-1);
    mintheta.push_back(-2);

    double dt=0.5;//discrete time interval dt
    int numInt = 10;//number of intervals/grind points on each parameter dimension
    int pre_kstep = 20;
    int test_set_start=1;
    int test_set_end = 20;
    for (int i=0;i<maxtheta.size();i++){
        double delta_temp = (maxtheta[i] - mintheta[i])/numInt/2.0;;
        delta_p.push_back(delta_temp);
    }
    // initialization:
    std::map<std::string, FnPtr> modelMap;
    string model="DubinsModel";
    modelMap[model] = DubinsModel_dynamic;


    
    vector <vector<double>> theta_list;
    vector <vector <int>> cell_table;
    
    cell_generator (delta_p.size(),numInt-1,cell_table);
    
    for (int i=0;i<cell_table.size();i++){
        vector <double> theta_centra;
        for (int j=0; j<cell_table[0].size(); j++){
            double temp=mintheta[j]+delta_p[j]*(2*cell_table[i][j]+1);
            theta_centra.push_back(temp);
        }
        theta_list.push_back(theta_centra);
    }
    
    //uniform dist
    vector <double> logPrior;
    double uniform_c=1.0;
    for (int i=0;i<mintheta.size();i++){
        uniform_c=uniform_c/numInt;
    }
    double loguniform_c=log(uniform_c);
    for (int i=0;i<cell_table.size();i++){
        logPrior.push_back(loguniform_c);
    }
    
    int sdim=4;//observed state variable
    
    ////read data
    vector <vector<double>> M,M_training,M_test;
    
    std::ifstream infile("../data/car_data.txt");
    std::string line;
    while (std::getline(infile, line))
    {
        vector <double> row;
        double a;
        std::istringstream iss(line);
        while((iss >> a)){// Read in from line stream
            row.push_back(a);
        }
        M.push_back(row);
    }

    //std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
    //std::cout.precision(2);
    //cout<<M[0][0]<<" "<<M[0][45]<<endl;
    //cout<<M.size()<<" "<< M[0].size()<<endl;
                 
//    int training_set_start=1;//test_set_start-50;
//    int training_set_end =50; //test_set_start;
//    int NT_train=training_set_end- training_set_start;
//
//    for(int i=training_set_start-1;i<training_set_end;i++){
//        M_training.push_back(M[i]);
//    }
    for(int i=test_set_start-1;i<test_set_end;i++){
        M_test.push_back(M[i]);
    }

    int idx_xx = 0;
    int idx_xy = 1;
    int idx_v = 2;
    int idx_psi = 3;
    
    
    // Calc cov matrix
    vector <vector<double>> cov(sdim,vector<double> (sdim)), inv_cov(sdim,vector<double> (sdim));
    //vector <vector<double>> error_train(sdim,vector<double> (NT_train-1));
    cov[0][0]=0.03;
    cov[0][1]=0.001;
    cov[1][0]=cov[0][1];
    cov[0][2]=0.01;
    cov[2][0]=cov[0][2];
    cov[2][2]=cov[0][2];
    cov[0][3]=0.01;
    cov[3][0]=cov[0][3];
    cov[1][1]=0.03;
    cov[1][2]=0.02;
    cov[2][1]=cov[1][2];
    cov[1][3]=0.01;
    cov[3][1]=cov[1][3];
    cov[2][2]=0.03;
    cov[2][3]=0.001;
    cov[3][2]=cov[2][3];
    cov[3][3]=0.03;
    
    double pi=3.14159265;

//    for(int i=0;i<NT_train-1;i++){
//        double uv_temp,upsi_temp;
//        vector<double> theta_temp, next_x(sdim),cur_x(sdim);
//        uv_temp = ( M_training[i+1][idx_v]- M_training[i][idx_v])/dt;
//        upsi_temp = ( M_training[i+1][idx_psi]- M_training[i][idx_psi])/dt;
//
//        theta_temp.push_back(uv_temp);
//        theta_temp.push_back(upsi_temp);
//
//        cur_x[0] = M_training[i][idx_xx];
//        cur_x[1] = M_training[i][idx_xy];
//        cur_x[2] = M_training[i][idx_v];
//        cur_x[3] = M_training[i][idx_psi];
//        next_x[0] = M_training[i+1][idx_xx];
//        next_x[1] = M_training[i+1][idx_xy];
//        next_x[2] = M_training[i+1][idx_v];
//        next_x[3] = M_training[i+1][idx_psi];
//
//        vector<double> simNext_x = modelMap[model](cur_x,theta_temp,dt);
//        //cout<<W_temp<<endl;
//        for(int j=0;j<sdim;j++) error_train[j][i]= next_x[j] - simNext_x[j];
//    }
//    //compute the covariance of error
//    CalculateVariane_matrix(error_train,cov);
    //compute the inverse of covariance
    inverse_cov(cov,inv_cov);
//    cout<<cov[0][0]<<" "<<cov[0][1]<<" "<<cov[0][2]<<" "<<cov[0][3]<<" "<<endl;
//    cout<<cov[1][0]<<" "<<cov[1][1]<<" "<<cov[1][2]<<" "<<cov[1][3]<<" "<<endl;
//    cout<<cov[2][0]<<" "<<cov[2][1]<<" "<<cov[2][2]<<" "<<cov[2][3]<<" "<<endl;
//    cout<<cov[3][0]<<" "<<cov[3][1]<<" "<<cov[3][2]<<" "<<cov[3][3]<<" "<<endl;

//    cout<<inv_cov[0][0]<<" "<<inv_cov[0][1]<<" "<<inv_cov[0][2]<<" "<<inv_cov[0][3]<<" "<<endl;
//    cout<<inv_cov[1][0]<<" "<<inv_cov[1][1]<<" "<<inv_cov[1][2]<<" "<<inv_cov[1][3]<<" "<<endl;
//    cout<<inv_cov[2][0]<<" "<<inv_cov[2][1]<<" "<<inv_cov[2][2]<<" "<<inv_cov[2][3]<<" "<<endl;
//    cout<<inv_cov[3][0]<<" "<<inv_cov[3][1]<<" "<<inv_cov[3][2]<<" "<<inv_cov[3][3]<<" "<<endl;
    int NT=test_set_end- test_set_start;
    vector <vector<double>> x(NT+1,vector<double> (sdim));
    for(int t=0;t<NT;t++){
        //data of state at t
        x[t][0] = M_test[t][idx_xx];
        x[t][1] = M_test[t][idx_xy];
        x[t][2] = M_test[t][idx_v];
        x[t][3] = M_test[t][idx_psi];
        //data of state at t+1
        x[t+1][0] = M_test[t+1][idx_xx];
        x[t+1][1] = M_test[t+1][idx_xy];
        x[t+1][2] = M_test[t+1][idx_v];
        x[t+1][3] = M_test[t+1][idx_psi];
    }


    //compute posterior p(omega(t)|x(t+1))
    int gridsize=logPrior.size();
    vector<vector<double>> logPosterior_list;
    for(int t=0;t<NT;t++){
        vector<double> logPosterior(gridsize);
        computePosterior(numInt,sdim,model,dt,t,x[t],x[t+1],inv_cov, theta_list, modelMap, logPrior, logPosterior);
        logPosterior_list.push_back(logPosterior);
        logPrior=logPosterior;
        //para.prior=(1-para.epsilon)*post+para.epsilon*para.prior_initial;
    }

//    for(int i=0;i<gridsize;i++){
//        cout<<logPosterior_list[NT-1][i]<<endl;
//    }

    std::ofstream output_logPosterior_list("logPosterior_list.txt");
    for (int i=0;i<logPosterior_list.size();i++)
    {
        for (int j=0;j<logPosterior_list[0].size();j++)
        {
            output_logPosterior_list << logPosterior_list[i][j] << " ";
        }
        output_logPosterior_list << "\n";
    }
    output_logPosterior_list.close();
    
    std::ofstream output_theta_list("theta_list.txt");
    for (int i=0;i<theta_list.size();i++)
    {
        for (int j=0;j<theta_list[0].size();j++)
        {
            output_theta_list << theta_list[i][j] << " ";
        }
        output_theta_list << "\n";
    }
    output_theta_list.close();
    


}
