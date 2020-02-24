//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_PROBABILITYQUERYEVALUATOR_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_PROBABILITYQUERYEVALUATOR_HH

#include "MultivariatePoly.hh"
#include "DistributionInfo.hh"

namespace PolynomialForms {

    class ProbabilityQueryEvaluator {
    protected:
        MultivariatePoly mp;
        std::map<int, DistributionInfoPtr > distributionInfo;
        std::map<int, MpfiWrapper> distributionRanges;
        std::vector<MultivariatePoly> splitComponents;
        std::vector<MpfiWrapper> componentExpectations;
        std::vector<MpfiWrapper> componentRanges;
        std::vector<MpfiWrapper> componentVariances;
        MpfiWrapper polynomialExpectation;
        MpfiWrapper expect;
        double t;
    public:
        ProbabilityQueryEvaluator(MultivariatePoly const & mp_, std::map<int, DistributionInfoPtr > const & dInfo_):
        mp(mp_), distributionInfo(dInfo_){
            expect = mp.expectation(distributionInfo);
            MpfiWrapper cTerm = mp.constantIntvl();
            mp.setConst(0);
            t = - cTerm.upper(); // if p + [l,u] >= 0  then p >= -u Pr(p + [l,u] >= 0) <= P(p >= -u) <= ..
            std::cout << "Debug: Converting p +  " << cTerm << " >= s into p >= s + " << t << std::endl;
            for (auto p: distributionInfo){
                distributionRanges.insert(std::make_pair(p.first, p.second->getRange()));
            }
        }

        void printPolyStats();
        void separatePolynomialIntoComponents();
        double computeChebyshevBounds() const;
        // Note this computes the log_e(bound). The user should compute exp(bound) returned by this function.
        double computeChernoffBound() const;
        // Note this computes the log_e(bound). The user should compute exp(bound) returned by this function.
        double computeBernsteinBound() const;
        double computeFourthMomentBound() const;
        void computeBestUpperTailBounds(double t);


    };

};

#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_PROBABILITYQUERYEVALUATOR_HH
