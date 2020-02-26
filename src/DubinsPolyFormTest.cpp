//
// Created by Sriram Sankaranarayanan on 2/25/20.
//

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include "DubinsPolyForm.hh"

bool debug = false;

void dubinsGenerateReachData(std::vector<double> u_vValues, std::vector<double> u_psiValues,
        double v0Lower, double v0Upper, double v0Delta, double deltaT, int numSteps, int maxDegree, ostream & outf){
    for (auto u_v: u_vValues){
        for (auto u_psi: u_psiValues){
            for (double v0 = v0Lower; v0 < v0Upper; v0 = v0 + v0Delta){
                std::cout << "u_v = " << u_v << ", u_psi = " << u_psi << ", v0 in " <<  v0 << ", " <<v0+v0Delta << std::endl;
                DubinsPolyForm dpf(u_v, u_psi, v0, v0+v0Delta, deltaT, numSteps, maxDegree);
                dpf.computeNStepForm(debug);
                dpf.printExpectationsAndRanges();
                dpf.printSpreadsheetRow(outf);
                //dpf.testWithSimulations(10000);

            }
        }
    }
}

int main(int argc, char * argv[]){
    int numSteps = 5;
    int maxDegree = 2;
    double deltaT = 0.5;
    std::string outFileName = "output.csv";
    if (argc > 1){
        outFileName = std::string(argv[1]);
    }
    if (argc > 2){
        numSteps = atoi(argv[2]);
        std::cout << "Steps: " << numSteps << std::endl;
    }
    if (argc > 3){
        maxDegree = atoi(argv[3]);
        std::cout << "Max Degree: " << maxDegree << std::endl;
    }

    ofstream outFileHandle(outFileName);
    std::vector<double> u_vValues = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9};
    std::vector<double> u_psiValues = {-1.8, -1.4, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.4, 1.8};
    double v0Lower = 10.0;
    double v0Upper = 20;
    double v0Delta = 0.5;
    auto start = chrono::high_resolution_clock::now();
    dubinsGenerateReachData(u_vValues, u_psiValues, v0Lower, v0Upper, v0Delta, deltaT,numSteps,maxDegree, outFileHandle);
    std::cout << "Output CSV in " << outFileName;
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    std::cout << "Time Taken: " << time_taken << std::endl;
    outFileHandle.close();
    return 1;
}