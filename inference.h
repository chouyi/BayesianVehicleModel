#include "dynamic.h"
# include "linpack_d.hpp"
#include <vector>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <unordered_map>
#include <map>
#include <iomanip>
//using namespace std;
typedef vector<double> (*FnPtr)(vector<double>,vector<double>,double);
#ifndef inference
#define inference
double CalculateMean(vector<double> & value);
double CalculateVariane(vector<double> & value);
double covariance(vector<double> & v1, vector<double> & v2);
void inverse_cov(vector<vector<double>> & cov,vector<vector<double>> & inv_cov);
void CalculateVariane_matrix(vector<vector<double>> & data,vector<vector<double>> & cov);
void getsubcell (const int m,const int k_order,vector<vector<int> > & v);
void cell_generator (const int n,const int kmax, vector<vector<int> > & table );
double log_mvnpdf(const vector<double> & mu,const vector<double> & x,const vector <vector<double>>& inv_cov);
void computePosterior(const int numInt,const int sdim,const string model,const double dt,const int t,const vector<double> & cur_x,const vector<double> & next_x,const vector <vector<double>>& inv_cov,const vector <vector<double>> theta_list,map<string, FnPtr> & modelMap,vector<double> & logPrior,vector<double> & logPosterior);
#endif
