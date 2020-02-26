//
// Created by Sriram Sankaranarayanan on 2/24/20.
//
#include <fstream>
#include "CTModelGenerateReachabilityData.hh"

void ctModelGenerateReachabilityData(std::vector<double> omegaValues,
                                     double vLow, double vHi, double vDelta,
                                     double deltaT,
                                     int numSteps,  int maxDegree, std::string outputFileName){
    // We need to iterate through the values of omega
    //   In turn, we need to iterate between vLow and vHi through steps of vDelta
    ofstream outFileHandle(outputFileName);
    for (double omega0: omegaValues){
        double vx0 = vLow;
        for (; vx0 < vHi; vx0 += vDelta){
            CTModelPolynomialForm ctmp(omega0, vx0, vx0+vDelta, -0.05, 0.05, deltaT, numSteps, maxDegree );
            std::cout << "omega: " << omega0 << std::endl;
            std::cout << "vx0 : " << vx0 << " , " << vx0 + vDelta << std::endl;
            ctmp.computeNStepForm();
            ctmp.printExpectationsAndRanges();
            outFileHandle << omega0 << "," << vx0 << "," << vx0 + vDelta ;
            ctmp.printSpreadsheetRow(outFileHandle);
            outFileHandle << std::endl;
            std::cout << "---------------------------------" << std::endl;
        }
    }
    outFileHandle.close();

}
