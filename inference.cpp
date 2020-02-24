#include "inference.h"

double CalculateMean(vector<double> & value)
{
    double sum = 0;
    int n=value.size();
    for(int i = 0; i < n; i++)
        sum += value[i];
    return (sum / n);
}



double CalculateVariane(vector<double> & value)
{
    double mean = CalculateMean(value);
    double temp = 0;
    int n= value.size();
    for(int i = 0; i < n; i++)
    {
        temp += pow((value[i] - mean),2) ;
        
    }
    return temp / n;
}

double covariance(vector<double> & v1, vector<double> & v2)
{
    int n=v1.size();
    int m=v2.size();
    if(n!=m) {cout<<"error:size difference"<<endl;}
    double sum = 0;
    for(int i = 0; i < n; i++)
        sum = sum + (v1[i] - CalculateMean(v1)) *
        (v2[i] - CalculateMean(v2));
    return sum / (n - 1);
}
void CalculateVariane_matrix(vector<vector<double>> & data,vector<vector<double>> & cov){
    //error_train(sdim,vector<double> (N)
    int d = data.size();
    for(int i=0;i<d;i++){
        for(int j=i; j<d;j++){
            cov[i][j]=covariance(data[i],data[j]);
            if(i!=j) cov[j][i]=cov[i][j];
        }
    }
}
void inverse_cov(vector<vector<double>> & cov,vector<vector<double>> & inv_cov){
    int LDA=cov.size();
    int N=cov[0].size();
    double *a = new double[LDA*N];
    for (int i = 1; i <= N; i++ )
    {
        for (int j = 1; j <= N; j++ ){
            a[i-1+(j-1)*LDA]=cov[i-1][j-1];
        }
    }

//    //   double a[] = {1, 3 ,5 ,2 ,5 , 6 , 3 , 5 , 9 };
    double det[2];
    int i;
    int info;
    int ipvt[N];
    int j;
    int job;
    double work[N];

    info = dgefa ( a, LDA, N, ipvt );

    if ( info != 0 )
    {
        cout << "  Error!  The matrix is nearly singular!\n";
        //return;
    }

    job = 11;
    dgedi ( a, LDA, N, ipvt, det, work, job );
    // << "  The determinant = " << det[0] << " * 10^"<< det[1] << "\n";
//    cout << "  The inverse matrix:\n";
//    cout << "\n";
//
//    for (int i = 1; i <= N; i++ )
//    {
//        for (int j = 1; j <= N; j++ )
//        {
//            cout << "  " << setw(12) << a[i-1+(j-1)*LDA];
//        }
//        cout << "\n";
//    }
    for (int i = 1; i <= N; i++ )
    {
        for (int j = 1; j <= N; j++ )
        {
            inv_cov[i-1][j-1] = a[i-1+(j-1)*LDA];
        }

    }

}


void getsubcell (const int m,const int k_order,vector<vector<int> > & v){
    if((m-1)>0){
        
        for (int i=0; i<=k_order;i++){
            
            vector<vector<int> > temp_v;
            getsubcell(m-1,k_order,temp_v);
            
            for(int j=0;j<temp_v.size();j++){
                temp_v[j].push_back(i);
            }
            for (int j=0;j<temp_v.size();j++){
                v.push_back(temp_v[j]);
            }
            
        }
    }
    else{
        for (int i=0; i<=k_order;i++){
            vector <int> v1;
            
            v1.push_back(i);
            v.push_back(v1);
        }
    }
    
}

void cell_generator (const int n,const int kmax, vector<vector<int> > & table ){
    //n: the number of variable
    //kmax: the maximum of interval's index(0,1,2,...,kmax)
    if(table.size() > 0){
        return;
    }
    
    vector<vector<int> > temp_table;
    getsubcell(n,kmax,temp_table);
    
    for (int j=0;j<temp_table.size();j++){
        vector<int> row;
        for (int k=temp_table[0].size();k>0;k--){
            row.push_back(temp_table[j][k-1]);
        }
        table.push_back(row);
    }
    temp_table.clear();
    
    
    
}

double log_mvnpdf(const vector<double> & mu,const  vector<double> & x,const vector <vector<double>>& inv_cov){
    //log Likelihood function : unnormalized value
    int n=mu.size();
    vector<double> dx(n);
    for(int i=0;i<n;i++){
        dx[i]=x[i]-mu[i];
    }
    
    vector<double> v(n,0.0);
    for(int i=0;i<inv_cov.size();i++){
        for(int j=0;j<n;j++){
            v[i]=v[i]+inv_cov[i][j]*dx[j];
        }
    }
//    //cout<<inv_cov[0][0]<<" "<<inv_cov[0][1]<<" "<<inv_cov[0][2]<<" "<<inv_cov[0][3]<<endl;
//    cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<" "<<v[3]<<endl;
//    cout<<v[0]<<" "<<dx[1]<<" "<<dx[2]<<" "<<dx[3]<<endl;
    double res=0;
    for(int i=0;i<n;i++){
        res+=dx[i]*v[i];
    }
    

    res=-res/2.0;
    return res;
}

void computePosterior(const int numInt,const int sdim,const string model,const double dt,const int t,const vector<double> & cur_x,const vector<double> & next_x,const vector <vector<double>>& inv_cov,const vector <vector<double>> theta_list,map<string, FnPtr> & modelMap,vector<double> & logPrior,vector<double> & logPosterior){
    
    for(int i=0;i<theta_list.size();i++){
        
        vector<double> simNext_x = modelMap[model](cur_x,theta_list[i],dt);
        double logp=0;
//        for(int j=0;j<sdim;j++){
//            logp=logp-pow((next_x[j]-simNext_x[j]),2)/(2*cov[j][j]);
//        }
        double pi=3.14159265;

        //cout<<exp(log_mvnpdf(next_x,simNext_x,inv_cov))/sqrt(pow(2*pi,4)*0.000108736)<<endl;
        //if(t==3){cout<<exp(log_mvnpdf(next_x,simNext_x,inv_cov))<<endl;}
        logp=logp+log_mvnpdf(simNext_x,next_x,inv_cov);
        logPosterior[i]= logp+logPrior[i];

        
    }
    
    
    
    //normalization
    double normC=0.0,lognormC;
    for(int i=0;i<logPosterior.size();i++){
        if (isfinite(logPosterior[i])== false){ logPosterior[i]=-1e+12;}
        normC+=exp(logPosterior[i]);
    }
    lognormC=log(normC);
    if (isfinite(lognormC)== false){ lognormC=-1e+12;}
    for(int i=0;i<logPosterior.size();i++){
        logPosterior[i]-=lognormC;
    }

    
    
}
