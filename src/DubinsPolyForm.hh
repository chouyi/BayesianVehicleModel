//
// Created by Sriram Sankaranarayanan on 2/25/20.
//

#ifndef BAYESIANVEHICLEMODEL_DUBINSPOLYFORM_HH
#define BAYESIANVEHICLEMODEL_DUBINSPOLYFORM_HH


#include <MpfiWrapper.hh>
#include <DistributionInfo.hh>
#include <MultivariatePoly.hh>
#include <ProbabilityQueryEvaluator.hh>

using namespace PolynomialForms;

class DubinsPolyForm{
protected:
    MpfiWrapper u_v0;
    MpfiWrapper u_psi0;
    MpfiWrapper v0Range;
    MpfiWrapper psi0Range;
    MpfiWrapper deltaT;

    int numSteps;
    int randomVarID;
    int maxDegree;


    MultivariatePoly x;
    MultivariatePoly y;
    MultivariatePoly v;
    MultivariatePoly psi;
    MultivariatePoly u_v;
    MultivariatePoly u_psi;

    std::map<int, DistributionInfoPtr > distrInfo;
    std::map<int, MpfiWrapper> env;

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


public:

    DubinsPolyForm(double uv0, double upsi0, double v0Lower, double v0Upper,
            double deltaT_, int numSteps_, int maxDegree_):
    u_v0(uv0),
    u_psi0(upsi0),
    v0Range(v0Lower, v0Upper),
    psi0Range(0.0),
    deltaT(deltaT_),
    numSteps(numSteps_),
    randomVarID(0),
    maxDegree(maxDegree_){};

    void computeNStepForm(bool debug=false);
    void printExpectationsAndRanges();
    void printSpreadsheetRow(ostream & outFileHandle);
    void testWithSimulations(int numSims);
    std::vector<double> runOneSimulation();



};
#endif //BAYESIANVEHICLEMODEL_DUBINSPOLYFORM_HH
