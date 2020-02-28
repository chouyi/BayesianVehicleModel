//
// Created by Sriram Sankaranarayanan on 2/24/20.
//

#include "CTModelPolyForm.hh"
#include "CTModelGenerateReachabilityData.hh"
#include <chrono>

int test1(){
    CTModelPolynomialForm ctmp(0.2, 12.1, 12.3, -0.1, 0.1, 0.4, 5);
    ctmp.computeNStepForm(true);
    double expectX = 0.0;
    double expectY = 0.0;
    double maxX = -1000000.0;
    double minX = 10000000.0;
    double maxY = -100000000.0;
    double minY = 100000000.0;
    int numSims = 10000;
    for (int j = 0; j < numSims; ++j){
        std::vector<double> xy = ctmp.simulateModel();
     //   std::cout << "x:" <<xy[0] << "y: " << xy[1] << std::endl;
        expectX = expectX + xy[0];
        expectY = expectY + xy[1];
        maxX = max(xy[0], maxX);
        minX = min(xy[0], minX);
        maxY = max(xy[1], maxY);
        minY = min(xy[1], minY);
    }

    ctmp.printExpectationsAndRanges();
    std::cout << "After " << numSims << " simulations: " << std::endl;
    std::cout << "E(x): " << expectX/(double) numSims << std::endl;
    std::cout << "E(y): " << expectY/(double) numSims << std::endl;
    std::cout << "Empirical Range x: [" << minX << "," << maxX << "]" << std::endl;
    std::cout << "Empirical Range y: [" << minY << "," << maxY << "]" << std::endl;
    return 1;
}

int main(int argc, char * argv[]){
    int numSteps = 5;
    int maxDegree = 1;
    string outFileName = "output.csv";
    double vx0Lower = 8.0;
    double vx0Upper = 22.0;
    double vx0Delta = 0.1;
    if (argc > 1){
        outFileName = string(argv[1]);
        std::cout << "Output File: " << outFileName << std::endl;
    }
    if (argc > 2){
        numSteps = atoi(argv[2]);
        std::cout << "Steps: " << numSteps << std::endl;
    }
    if (argc > 3){
        maxDegree = atoi(argv[3]);
        std::cout << "Max Degree: " << maxDegree << std::endl;
    }

    if (argc > 5){
        vx0Lower = atof(argv[4]);
        vx0Upper = atof(argv[5]);
        std::cout << "Velocity range: " << vx0Lower << "," << vx0Upper << std::endl;
    }

    std::vector<double> omega0Values = { -0.171, -0.153,-0.135, -0.117, -0.099,   -0.081,  -0.063, -0.045,
                                         -0.027,  -0.009,0.009, 0.027, 0.045, .063, 0.081,  0.099, 0.117,
                                         0.135,  0.153, 0.171 };

    auto start = chrono::high_resolution_clock::now();
    ctModelGenerateReachabilityData(omega0Values, vx0Lower, vx0Upper, vx0Delta, 0.4, numSteps, maxDegree, outFileName);
    std::cout << "Output CSV in " << outFileName;
    auto end = chrono::high_resolution_clock::now();
    double time_taken =
            chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;
    std::cout << "Time Taken: " << time_taken << std::endl;
    return 1;
}