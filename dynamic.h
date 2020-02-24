#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

#ifndef dynamic
#define dynamic
vector<double> CTmodel_dynamic(vector<double> cur_x,vector<double> theta,double dt);
vector<double> DubinsModel_dynamic(vector<double> x,vector<double> theta,double dt);
#endif
