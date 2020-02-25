
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

const string ctModelDataFile = "../data/UAV_data.txt";

//template<class T>
//void display(T A[N][N])
//{
//    for (int i=0; i<N; i++)
//    {
//        for (int j=0; j<N; j++)
//            cout << A[i][j] << " ";
//        cout << endl;
//    }
//}
int main(){
    vector<double> maxtheta,mintheta,delta_p;
    maxtheta.push_back(0.18);
    mintheta.push_back(-0.18);

    double dt=0.4;//discrete time interval dt
    int numInt = 20;//number of intervals/grind points on each parameter dimension
    int pre_kstep = 20;
    int test_set_start=4200;
    int test_set_end = 4300;
    for (int i=0;i<maxtheta.size();i++){
        double delta_temp = (maxtheta[i] - mintheta[i])/numInt/2.0;;
        delta_p.push_back(delta_temp);
    }
    // initialization:
    std::map<std::string, FnPtr> modelMap;
    string model="CTmodel";
    modelMap[model] = CTmodel_dynamic;


    
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
    std::ifstream infile(ctModelDataFile);
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
                 
    int training_set_start=test_set_start-1000;
    int training_set_end = test_set_start;
    int NT_train=training_set_end- training_set_start;

    for(int i=training_set_start-1;i<training_set_end;i++){
        M_training.push_back(M[i]);
    }
    for(int i=test_set_start-1;i<test_set_end;i++){
        M_test.push_back(M[i]);
    }

    int idx_xy = 3-1;
    int idx_xx = 4-1;
    int idx_vy = 6-1;
    int idx_vx = 7-1;
    
    
    // Calc cov matrix
    vector <vector<double>> cov(sdim,vector<double> (sdim)), error_train(sdim,vector<double> (NT_train-1)),inv_cov(sdim,vector<double> (sdim));
//    cov[0][0]=pow(1.188274283,2);
//    cov[1][1]=pow(1.174298356,2);
//    cov[2][2]=pow(1.996316081,2);
//    cov[3][3]=pow(1.128651877,2);
    double pi=3.14159265;

    for(int i=0;i<NT_train-1;i++){
        double W_temp;
        vector<double> theta_temp, next_x(sdim),cur_x(sdim);
        W_temp = atan2(M_training[i+1][idx_vy],M_training[i+1][idx_vx]) - atan2(M_training[i][idx_vy]  ,M_training[i][idx_vx]);
        if(W_temp>pi){
            W_temp = W_temp - ceil(W_temp/pi)*pi;
        }else if (W_temp < -pi){
            W_temp = W_temp + ceil(W_temp/pi)*pi;
        }
        theta_temp.push_back(W_temp);
        cur_x[0] = M_training[i][idx_xx];
        cur_x[1] = M_training[i][idx_vx];
        cur_x[2] = M_training[i][idx_xy];
        cur_x[3] = M_training[i][idx_vy];
        next_x[0] = M_training[i+1][idx_xx];
        next_x[1] = M_training[i+1][idx_vx];
        next_x[2] = M_training[i+1][idx_xy];
        next_x[3] = M_training[i+1][idx_vy];

        vector<double> simNext_x = CTmodel_dynamic(cur_x,theta_temp,dt);
        //cout<<W_temp<<endl;
        for(int j=0;j<sdim;j++) error_train[j][i]= next_x[j] - simNext_x[j];
    }
    //compute the covariance of error
    CalculateVariane_matrix(error_train,cov);
    //compute the inverse of covariance
    inverse_cov(cov,inv_cov);
//    cout<<inv_cov[0][0]<<" "<<inv_cov[0][1]<<" "<<inv_cov[0][2]<<" "<<inv_cov[0][3]<<" "<<endl;
//    cout<<inv_cov[1][0]<<" "<<inv_cov[1][1]<<" "<<inv_cov[1][2]<<" "<<inv_cov[1][3]<<" "<<endl;
//    cout<<inv_cov[2][0]<<" "<<inv_cov[2][1]<<" "<<inv_cov[2][2]<<" "<<inv_cov[2][3]<<" "<<endl;
//    cout<<inv_cov[3][0]<<" "<<inv_cov[3][1]<<" "<<inv_cov[3][2]<<" "<<inv_cov[3][3]<<" "<<endl;
    int NT=test_set_end- test_set_start;
    vector <vector<double>> x(NT+1,vector<double> (sdim));
    for(int t=0;t<NT;t++){
        //data of state at t
        x[t][0] = M_test[t][idx_xx];
        x[t][1] = M_test[t][idx_vx];
        x[t][2] = M_test[t][idx_xy];
        x[t][3] = M_test[t][idx_vy];
        //data of state at t+1
        x[t+1][0] = M_test[t+1][idx_xx];
        x[t+1][1] = M_test[t+1][idx_vx];
        x[t+1][2] = M_test[t+1][idx_xy];
        x[t+1][3] = M_test[t+1][idx_vy];
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

    
    



}
