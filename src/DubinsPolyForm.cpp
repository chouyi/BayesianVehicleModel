//
// Created by Sriram Sankaranarayanan on 2/25/20.
//

#include "DubinsPolyForm.hh"
#include <cmath>

void DubinsPolyForm::initializeForms() {
    x.setConst(0.0);
    y.setConst(0.0);
    v.setConst(v0Range);
    psi.setConst(psi0Range);
    u_v.setConst(u_v0);
    u_psi.setConst(u_psi0);
}

void DubinsPolyForm::computeOneStep() {
    /*
     * next_x[0]=x +h* v*cos(psi)+(h*h/2)*u_v*cos(psi);
    next_x[1]=x[1]+h*x[2]*sin(x[3])+(h*h/2)*u_v*sin(x[3]);
    next_x[2]=x[2]+h*u_v;
    next_x[3]=x[3]+h*u_psi;
     */
    MultivariatePoly cosPsi = psi.cosine(env);
    MultivariatePoly sinPsi = psi.sine(env);

    MultivariatePoly vCosPsi = v.multiply(cosPsi);
    MultivariatePoly vSinPsi = v.multiply(sinPsi);
    MultivariatePoly u_vCos = u_v.multiply(cosPsi);
    MultivariatePoly u_vSin = u_v.multiply(sinPsi);

    x.scaleAndAddAssign(deltaT, vCosPsi);
    x.scaleAndAddAssign(deltaT* deltaT * MpfiWrapper(0.5), u_vCos );
    x.centerAssign(env);
    x.truncateAssign(maxDegree, env);

    y.scaleAndAddAssign(deltaT, vSinPsi);
    y.scaleAndAddAssign(deltaT * deltaT * MpfiWrapper(0.5), u_vSin);
    y.centerAssign(env);
    y.truncateAssign(maxDegree, env);

    v.scaleAndAddAssign(deltaT, u_v);
    psi.scaleAndAddAssign(deltaT, u_psi);

    int u_vChange = createUniform(-0.25, 0.25);
    MultivariatePoly u_vChangePoly(1.0, u_vChange);
    u_v.scaleAndAddAssign(1.0, u_vChangePoly);

    int u_psiChange = createUniform(-0.02, 0.02);
    MultivariatePoly u_psiChangePoly(1.0, u_psiChange);
    u_psi.scaleAndAddAssign(1.0, u_psiChangePoly);
}




void DubinsPolyForm::computeNStepForm(bool debug){
    std::map<int, string> dummy_printer;
    initializeForms();
    for (int i=0; i < numSteps; ++i){
        if (debug){
            std::cout << "Step # " << i << std::endl;
        }
        computeOneStep();
        if (i % 5 == 0){
            MpfiWrapper psiRng = psi.evaluate(env);
            psiRng= intersect(psiRng, MpfiWrapper(-3.1415, 3.1415));
            psi = MultivariatePoly(psiRng);

        }

        if (debug){
            std::cout << "x = " << std::endl;
            x.prettyPrint(std::cout, dummy_printer);

            std::cout << "y = " << std::endl;
            y.prettyPrint(std::cout, dummy_printer);

        }
    }
}

void DubinsPolyForm::printExpectationsAndRanges() {
    MpfiWrapper ex = x.expectation(distrInfo);
    MpfiWrapper ey  = y.expectation(distrInfo);
    MpfiWrapper rngx = x.evaluate(env);
    MpfiWrapper rngy = y.evaluate(env);
    std::cout << rngx.lower() << "<= x <= " << rngx.upper() << std::endl;
    std::cout << rngy.lower() << "<= y <= " << rngy.upper() << std::endl;
    MultivariatePoly diff(x);
    diff.scaleAndAddAssign(-1.0, y);
    MpfiWrapper rng_x_minus_y = diff.evaluate(env);
    std::cout << rng_x_minus_y.lower() << "<= x - y <= " << rng_x_minus_y.upper() << std::endl;
    MultivariatePoly sumPoly(x);
    sumPoly.scaleAndAddAssign(1.0, y);
    MpfiWrapper rng_x_plus_y = sumPoly.evaluate(env);
    std::cout << rng_x_plus_y.lower() << "<= x + y <= " << rng_x_plus_y.upper() << std::endl;
}

void DubinsPolyForm::printSpreadsheetRow(ostream & outFileHandle) {
    MpfiWrapper rngx = x.evaluate(env);
    MpfiWrapper rngy = y.evaluate(env);
    MultivariatePoly diff(x);
    diff.scaleAndAddAssign(-1.0, y);
    MpfiWrapper rng_x_minus_y = diff.evaluate(env);
    rng_x_minus_y = intersect(rng_x_minus_y, rngx-rngy);
    MultivariatePoly sumPoly(x);
    sumPoly.scaleAndAddAssign(1.0, y);
    MpfiWrapper rng_x_plus_y = sumPoly.evaluate(env);
    rng_x_plus_y = intersect(rng_x_plus_y, rngx+rngy);
    outFileHandle << "," << rngx.lower() << "," << rngx.upper();
    outFileHandle << "," << rngy.lower() << "," << rngy.upper();
    outFileHandle << "," << rng_x_minus_y.lower() << "," << rng_x_minus_y.upper();
    outFileHandle << "," << rng_x_plus_y.lower() << "," << rng_x_plus_y.upper();
}

std::vector<double> DubinsPolyForm::runOneSimulation() {
    double x = 0;
    double y = 0;
    double psi = 0;
    double v =v0Range.getRandomSample();
    double u_v = median(u_v0);
    double u_psi = median(u_psi0);
    double dt = median(deltaT);
    MpfiWrapper u_vChange(-0.25, 0.25);
    MpfiWrapper u_psiChange(-0.02, 0.02);
    for (int i =0; i < numSteps; ++i){
        x  = x + dt * v * cos(psi) + 0.5 * dt * dt * u_v * cos(psi);
        y = y + dt * v * sin(psi) + 0.5 * dt * dt * u_v * sin(psi);
        v = v + dt * u_v;
        psi = psi + dt * u_psi;
        u_v = u_v + u_vChange.getRandomSample();
        u_psi = u_psi + u_psiChange.getRandomSample();
    }
    return std::vector<double>({x,y});
}


void DubinsPolyForm::testWithSimulations(int numSims){
    double xMin = 1000000.0;
    double xMax = -xMin;
    double yMin = 10000000.0;
    double yMax = -yMin;

    double xMinusYMin = xMin;
    double xMinusYMax = xMax;

    double xPlusYMin = xMin;
    double xPlusYMax = xMax;


    for (int j = 0; j < numSims; ++j){
        auto xy = runOneSimulation();
        xMin = min(xMin, xy[0]);
        xMax = max(xMax, xy[0]);

        yMin = min(yMin, xy[1]);
        yMax = max(yMax, xy[1]);

        xMinusYMin = min(xMinusYMin, xy[0] - xy[1]);
        xMinusYMax = max(xMinusYMax, xy[0] - xy[1]);

        xPlusYMin = min(xPlusYMin, xy[1] + xy[0]);
        xPlusYMax= max(xPlusYMax, xy[1] + xy[0]);
    }
    std::cout << "After " << numSims << " simulations: " << std::endl;
    std::cout << "x in " << xMin << "," << xMax << std::endl;
    std::cout << "y in " << yMin << "," << yMax << std::endl;
    std::cout << "x - yin " << xMinusYMin << "," << xMinusYMax << std::endl;
    std::cout << "x + y in " << xPlusYMin << "," << xPlusYMax << std::endl;

}
