//
// Created by Sriram Sankaranarayanan on 2/13/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_MODELPARSER_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_MODELPARSER_HH

#include "DistributionInfo.hh"
#include "SystemDescription.hh"
#include "ExprAST.hh"

using namespace PolynomialForms;
extern StochasticSystem * globalSystem;
extern int lineNum;

void parserMain(const char * fileName);
StateAbstractionPtr computeNSteps(int n, int maxDegree);



#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_MODELPARSER_HH
