//
// Created by Sriram Sankaranarayanan on 2/24/20.
//

#ifndef BAYESIANVEHICLEMODEL_CTMODELPOLYFORM_HH
#define BAYESIANVEHICLEMODEL_CTMODELPOLYFORM_HH

#include <MpfiWrapper.hh>
#include <DistributionInfo.hh>
#include <MultivariatePoly.hh>
#include <ProbabilityQueryEvaluator.hh>

using namespace PolynomialForms;

class CTModelPolynomialForm {
    protected:
    MpfiWrapper omega0;
    MpfiWrapper vx0;
    MpfiWrapper vy0;
    MpfiWrapper x0;
    MpfiWrapper y0;
    MpfiWrapper deltaT;


    std::map<int, DistributionInfoPtr> distrInfo;
    std::map<int, MpfiWrapper> env;
    MultivariatePoly x;
    MultivariatePoly y;
    MultivariatePoly vx;
    MultivariatePoly vy;
    MultivariatePoly omega;
    int numSteps;
    int randomVarID;

    void initializeForms();
    void computeOneStep();
    int createUniform(double a, double b){
        int varID = randomVarID;
        DistributionInfoPtr d1 = std::make_shared<UniformDistributionInfo>(MpfiWrapper(a,b));
        distrInfo.insert(std::make_pair(varID, d1));
        env.insert(std::make_pair(varID, MpfiWrapper(a,b)));
        randomVarID ++;
        return varID;
    }

    int createTruncGaussian(double mean, double sd, double a, double b){
        int varID = randomVarID;
        DistributionInfoPtr d1 = std::make_shared<TruncNormalDistributionInfo>
                (MpfiWrapper(a, b), mean, sd);
        distrInfo.insert(std::make_pair(varID, d1));
        env.insert(std::make_pair(varID, MpfiWrapper(a,b)));
        randomVarID ++;
        return varID;

    }

    void updateState( std::vector<double> & curState);

    public:

    CTModelPolynomialForm(double omega0_,
                        double vx0Lower_, double vx0Upper_,
                        double vy0Lower_, double vy0Upper_,
                        double deltaT_,
                        int numSteps_):
    omega0(omega0_),
    vx0(vx0Lower_, vx0Upper_),
    vy0(vy0Lower_, vy0Upper_),
    x0(0.0),
    y0(0.0),
    deltaT(deltaT_),
    numSteps(numSteps_),
    randomVarID(0)
    {}

    void computeNStepForm(bool debug=false);

    double boundProbability(double obstacleXLow, double obstacleXHi, double obstacleYLo, double obstacleYHi);

    void printExpectationsAndRanges();
    std::vector<double> simulateModel();

};


#endif //BAYESIANVEHICLEMODEL_CTMODELPOLYFORM_HH
