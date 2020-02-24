//
// Created by Sriram Sankaranarayanan on 2/24/20.
//
#include "ctModelGenerateReachabilityData.hh"

void ctModelGenerateReachabilityData(std::vector<double> omegaValues,
                                     double vLow, double vHi, double vDelta,
                                     double deltaT,
                                     int numSteps){
    // We need to iterate through the values of omega
    //   In turn, we need to iterate between vLow and vHi through steps of vDelta
    for (double omega0: omegaValues){
        double vx0 = vLow;
        for (; vx0 < vHi; vx0 += vDelta){
            CTModelPolynomialForm ctmp(omega0, vx0, vx0+vDelta, -0.05, 0.05, deltaT, numSteps );
            std::cout << "omega: " << omega0 << std::endl;
            std::cout << "vx0 : " << vx0 << " , " << vx0 + vDelta << std::endl;
            ctmp.computeNStepForm();
            ctmp.printExpectationsAndRanges();
            std::cout << "---------------------------------" << std::endl;
        }
    }

}
