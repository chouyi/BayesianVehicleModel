//
// Created by Sriram Sankaranarayanan on 2/24/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_POLYNOMIALFUNCTIONS_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_POLYNOMIALFUNCTIONS_HH
#include "MultivariatePoly.hh"

namespace PolynomialForms{
    /*-- Special functions for the coordinated turn models --*/
    MultivariatePoly computeSinXByX(MultivariatePoly const & p, std::map<int, MpfiWrapper> const &var_env,
            int resultDegree );
    MultivariatePoly computeCosXMinusOneByX(MultivariatePoly const & p, std::map<int, MpfiWrapper> const &var_env,
            int resultDegree );
};

#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_POLYNOMIALFUNCTIONS_HH
