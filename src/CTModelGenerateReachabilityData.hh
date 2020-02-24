//
// Created by Sriram Sankaranarayanan on 2/24/20.
//

#ifndef BAYESIANVEHICLEMODEL_CTMODELGENERATEREACHABILITYDATA_HH
#define BAYESIANVEHICLEMODEL_CTMODELGENERATEREACHABILITYDATA_HH
#include "CTModelPolyForm.hh"

void ctModelGenerateReachabilityData(std::vector<double> omegaValues,
                                        double vLow, double vHi,
                                        double vDelta,
                                        double deltaT,
                                        int numSteps);
#endif //BAYESIANVEHICLEMODEL_CTMODELGENERATEREACHABILITYDATA_HH
