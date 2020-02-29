//
// Created by Sriram Sankaranarayanan on 2/24/20.
//
#include <cmath>
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

MultivariatePoly specialSine(MultivariatePoly const & p, std::map<int, MpfiWrapper> const &var_env, int resultDegree) {
    /*--
     * Compute power series of sine around 0.0
     *   Just make sure that the power of x is one short
     *   Note that with remainder sin(x) =
     *       x - x^3/6  - cos(rng) x^4/24
     *    sin(x)/x = 1 - x^2/6 - cos(rng) x^3/24
     */

    MultivariatePoly res(1.0);
    MultivariatePoly tmp = p.squarePoly();
    MpfiWrapper rng = p.evaluate(var_env);
    MpfiWrapper cos_rng = cos(rng);

    //res = res -tmp/6
    res.scaleAndAddAssign(MpfiWrapper(-1.0/6.0), tmp);
    //tmp3 = p^3
    MultivariatePoly tmp3 = tmp.multiply(p);
    // res = res -1/24 cos(rng)*x^3
    res.scaleAndAddAssign(cos_rng * MpfiWrapper(-1.0/24), tmp3);
    res.centerAssign(var_env);
    res.truncateAssign(resultDegree, var_env);
    return res;
}

MultivariatePoly
specialCosine(MultivariatePoly const &p, std::map<int, MpfiWrapper> const &var_env, int resultDegree) {

    MultivariatePoly res(0.0);
    MpfiWrapper rng = p.evaluate(var_env);
    MpfiWrapper sin_rng = sin(rng);

    res.scaleAndAddAssign(MpfiWrapper(-0.5), p);
    MultivariatePoly tmp2 = p.squarePoly();
    MultivariatePoly tmp3 = p.multiply(tmp2);
    res.scaleAndAddAssign(MpfiWrapper(1.0/24), tmp3);
    MultivariatePoly tmp4 = tmp2.squarePoly();
    res.scaleAndAddAssign(MpfiWrapper(-1.0/120) * sin_rng, tmp4);
    res.centerAssign(var_env);
    res.truncateAssign(resultDegree, var_env);
    return res;
}


void CTModelPolynomialForm::computeOneStep(){
    MultivariatePoly omegaT = omega;
    omegaT.scaleAssign(deltaT); // omega * T

    MultivariatePoly sineOmegaTByOmegaT = specialSine(omegaT, env, maxDegree);
    sineOmegaTByOmegaT.scaleAssign(deltaT);

    MultivariatePoly cosOmegaTMinusOneByOmegaT = specialCosine(omegaT, env, maxDegree);
    cosOmegaTMinusOneByOmegaT.scaleAssign(deltaT);

    /*std::cout << "DEBUG: omegaT = " << std::endl;
    omegaT.prettyPrint(std::cout, std::map<int, string>());
    std::cout << std::endl;

    std::cout << "DEBUG: Sin(omega t)/omega = ";
    sineOmegaTByOmegaT.prettyPrint(std::cout, std::map<int, string>());
    std::cout << std::endl;
    std::cout << "DEBUG: (cos(omegat) - 1)/omega = ";
    cosOmegaTMinusOneByOmegaT.prettyPrint(std::cout, std::map<int, string>());
    std::cout << std::endl; */


    MultivariatePoly sineOmegaT = omegaT.sine(env);
    MultivariatePoly cosOmegaT = omegaT.cosine(env);

    /*-- Update for x --*/
    MultivariatePoly sinTimesVx = sineOmegaTByOmegaT.multiply(vx);
    //sinTimesVx.truncateAssign(maxDegree, env);
    MultivariatePoly cosTimesVy = cosOmegaTMinusOneByOmegaT.multiply(vy);
   // cosTimesVy.truncateAssign(maxDegree, env);
    x.scaleAndAddAssign(MpfiWrapper(1.0), sinTimesVx);
    x.scaleAndAddAssign(MpfiWrapper(1.0), cosTimesVy);
    x.truncateAssign(maxDegree,env);
    x.centerAssign(env);

    /*-- Update for y --*/
    MultivariatePoly cosTimesVx = cosOmegaTMinusOneByOmegaT.multiply(vx);
    //cosTimesVx.truncateAssign(maxDegree, env);
    MultivariatePoly sinTimesVy = sineOmegaTByOmegaT.multiply(vy);
   // sinTimesVy.truncateAssign(maxDegree, env);
    y.scaleAndAddAssign(-1.0, cosTimesVx);
    y.scaleAndAddAssign(1.0, sinTimesVy);
    y.truncateAssign(maxDegree, env);
    y.centerAssign(env);

    /*-- Update for vx --*/
    MultivariatePoly v11 = cosOmegaT.multiply(vx);
    //v11.truncateAssign(maxDegree, env);
    MultivariatePoly v12 = sineOmegaT.multiply(vy);
   // v12.truncateAssign(maxDegree, env);
    int vx_noiseID = createUniform(-0.2, 0.2);
    MultivariatePoly vx_noise(1.0,  vx_noiseID);
    vx = v11;
    vx.scaleAndAddAssign(-1.0, v12);
    vx.scaleAndAddAssign(1.0, vx_noise);
    vx.truncateAssign(maxDegree, env);
    vx.centerAssign(env);

    /*-- Update for vy --*/
    MultivariatePoly v21 = sineOmegaT.multiply(vx);
    MultivariatePoly v22 = cosOmegaT.multiply(vy);
    int vy_noiseID = createUniform(-0.2, 0.2);
    MultivariatePoly vy_noise(1.0,  vy_noiseID);
    vy = v21;
    vy.scaleAndAddAssign(1.0, v22);
    vy.scaleAndAddAssign(1.0, vy_noise);
    vy.truncateAssign(maxDegree, env);
    vy.centerAssign(env);

    /*-- Update for omega --*/
    int noiseTermID = createUniform(-0.045, 0.045);
    MultivariatePoly noiseTerm(1.0,  noiseTermID);
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
        if (i > 0 && i % 5 == 0){
            MpfiWrapper omegaRng = omega.evaluate(env);
            omegaRng= intersect(omegaRng, MpfiWrapper(omegaLow, omegaHi));
            omega = MultivariatePoly(omegaRng);

            MpfiWrapper vxRng = vx.evaluate(env);
            vx = MultivariatePoly(vxRng);

            MpfiWrapper vyRng = vx.evaluate(env);
            vy = MultivariatePoly(vyRng);


        }
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
        std::cout << "Range(x) = " << x.evaluate(env) << std::endl;
        std::cout << "Range(y) = " << y.evaluate(env) << std::endl;
        std::cout << "E(y) = " << ey << std::endl;


        // We will use some template

        return 1.0;
}

double sinxbyx(double x){
    if (fabs(x) <= 1E-02){
        return 1.0 - x*x/6.0;
    }
    return sin(x)/x;
}

double cosxminusonebyx(double x){
    if (fabs(x) <= 1E-02){
        return -0.5 * x + x*x*x/24.0;
    }
    return (cos(x)-1.0)/x;
}

void CTModelPolynomialForm::updateState( std::vector<double> & curState){
    double x = curState[0];
    double y = curState[1];
    double vx = curState[2];
    double vy = curState[3];
    double omega = curState[4];
    double dt = median(deltaT);

    curState[0] = x + dt* vx * sinxbyx(dt*omega) + dt * vy * cosxminusonebyx(dt * omega) ;
    curState[1] = y - dt * vx * cosxminusonebyx(dt * omega) + dt * vy * sinxbyx(dt * omega);
    curState[2] = cos(dt * omega) * vx - sin(dt*omega) * vy;
    curState[3] = sin(dt*omega)* vx + cos(dt*omega) * vy;
    curState[4] = curState[4] + (0.09 * (double) rand()/ (double) RAND_MAX -0.045);

}

std::vector<double> CTModelPolynomialForm::simulateModel() {
    std::vector<double> res;
    std::vector<double> curState;
    curState.push_back(x0.getRandomSample());
    curState.push_back(y0.getRandomSample());
    curState.push_back(vx0.getRandomSample());
    curState.push_back(vy0.getRandomSample());
    curState.push_back(omega0.getRandomSample());

    for (int i =0; i < numSteps; ++i){
        updateState(curState);
    }
    return vector<double>({curState[0], curState[1]});

}

void CTModelPolynomialForm::printExpectationsAndRanges() {
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

void CTModelPolynomialForm::printSpreadsheetRow(ostream & outFileHandle) {
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



