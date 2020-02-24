//
// Created by Sriram Sankaranarayanan on 2/24/20.
//

#include "CTModelPolyForm.hh"
#include "PolynomialFunctions.hh"

void CTModelPolynomialForm::initializeForms() {
    /*-- make up the initial distributions --*/
    // Let us initialize x, y : currently to (0,0)
    x.setConst(x0);
    y.setConst(y0);
    // v should be initialized to the interval vx0 and vy0 resp.
    vx.setConst(vx0);
    vy.setConst(vy0);
    // Initialize omega to omega0
    omega.setConst(omega0);
}



void CTModelPolynomialForm::computeOneStep(){
    MultivariatePoly omegaT = omega;
    omegaT.scaleAssign(deltaT); // omega * T
    MultivariatePoly sineOmegaTByOmegaT = computeSinXByX(omegaT, env, 1);
    sineOmegaTByOmegaT.scaleAssign(deltaT);

    MultivariatePoly cosOmegaTMinusOneByOmegaT = computeCosXMinusOneByX(omegaT, env, 1);
    cosOmegaTMinusOneByOmegaT.scaleAssign(deltaT);

    MultivariatePoly sineOmegaT = omegaT.sine(env);
    MultivariatePoly cosOmegaT = omegaT.cosine(env);

    /*-- Update for x --*/
    MultivariatePoly sinTimesVx = sineOmegaTByOmegaT.multiply(vx);
    sinTimesVx.truncateAssign(1, env);

    MultivariatePoly cosTimesVy = cosOmegaTMinusOneByOmegaT.multiply(vy);
    cosTimesVy.truncateAssign(1, env);

    x.scaleAndAddAssign(MpfiWrapper(1.0), sinTimesVx);
    x.scaleAndAddAssign(MpfiWrapper(1.0), cosTimesVy);
    x.centerAssign(env);

    /*-- Update for y --*/
    MultivariatePoly cosTimesVx = cosOmegaTMinusOneByOmegaT.multiply(vx);
    cosTimesVx.truncateAssign(1, env);
    MultivariatePoly sinTimesVy = sineOmegaTByOmegaT.multiply(vy);
    sinTimesVy.truncateAssign(1, env);
    y.scaleAndAddAssign(1.0, cosTimesVx);
    y.scaleAndAddAssign(1.0, sinTimesVy);
    y.centerAssign(env);

    /*-- Update for vx --*/
    MultivariatePoly v11 = cosOmegaT.multiply(vx);
    MultivariatePoly v12 = sineOmegaT.multiply(vy);
    v11.truncateAssign(1, env);
    v12.truncateAssign(1, env);
    vx = v11;
    vx.scaleAndAddAssign(-1.0, v12);
    vy.centerAssign(env);

    /*-- Update for vy --*/
    MultivariatePoly v21 = sineOmegaT.multiply(vx);
    MultivariatePoly v22 = cosOmegaT.multiply(vy);
    v21.truncateAssign(1, env);
    v22.truncateAssign(2, env);
    vy = v21;
    vy.scaleAndAddAssign(1.0, v22);
    vy.centerAssign(env);

    /*-- Update for omega --*/
    int noiseTermID = createUniform(-0.05, 0.05);
    MultivariatePoly noiseTerm(1.0, noiseTermID);
    omega.scaleAndAddAssign(1.0, noiseTerm);
}

void CTModelPolynomialForm::computeNStepForm(bool debug){
    std::map<int, string> dummy_printer;
    initializeForms();
    for (int i=0; i < numSteps; ++i){
        if (debug){
            std::cout << "Step # " << i << std::endl;
        }
        computeOneStep();
        if (debug){
            std::cout << "x = " << std::endl;
            x.prettyPrint(std::cout, dummy_printer);

            std::cout << "y = " << std::endl;
            y.prettyPrint(std::cout, dummy_printer);

        }
    }
}

double CTModelPolynomialForm::boundProbability(double obstacleXLow, double obstacleXHi, double obstacleYLo,
                                               double obstacleYHi) {
        // Assume we have called compute N Step Form before calling this function.
        MpfiWrapper ex = x.expectation(distrInfo);
        MpfiWrapper ey  = y.expectation(distrInfo);
        std::cout << "After " << numSteps << "steps " << std::endl;
        std::cout << "E(x) = " << ex << std::endl;
        std::cout << "E(y) = " << ey << std::endl;
        return 1.0;
}



