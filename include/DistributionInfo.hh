//
// Created by Sriram Sankaranarayanan on 2/10/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_DISTRIBUTIONINFO_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_DISTRIBUTIONINFO_HH

#include <string>
#include <map>
#include <vector>
#include "MpfiWrapper.hh"
#include "truncatedNormalCode.hh"

namespace PolynomialForms {

    class VariableDistributionInfo {
    protected:

        MpfiWrapper range;
    public:
        VariableDistributionInfo( MpfiWrapper range_):  range(range_){};
        virtual MpfiWrapper getRange() const { return range; };
        virtual MpfiWrapper getExpectation() const = 0;
        virtual MpfiWrapper getMoment(int j) const = 0;
        virtual ~VariableDistributionInfo() = default;
    };
    typedef std::shared_ptr<VariableDistributionInfo> DistributionInfoPtr;


    class TruncNormalDistributionInfo: public VariableDistributionInfo {
    protected:
        MpfiWrapper mean;
        MpfiWrapper sdev;
    public:
        TruncNormalDistributionInfo( MpfiWrapper range_, MpfiWrapper mean_, MpfiWrapper sdev_):
        VariableDistributionInfo(range_),
        mean(mean_),
        sdev(sdev_) {};

        virtual MpfiWrapper getExpectation() const { return mean; };

        virtual MpfiWrapper getMoment(int j) const {
            if (j == 1) { return getExpectation(); }
            double a = range.lower();
            double b = range.upper();
            double mu = median(mean);
            double sigma = median(sdev);
            return MpfiWrapper(truncated_normal_ab_moment (j, mu, sigma, a, b));
        }

    };


    class UniformDistributionInfo: public VariableDistributionInfo {
    protected:
        MpfiWrapper mean;
    public:
        UniformDistributionInfo(MpfiWrapper range_):
                VariableDistributionInfo(range_),
                mean((range_.upper() + range_.lower())/2.0)
                {};

        virtual MpfiWrapper getExpectation() const { return mean; };

        virtual MpfiWrapper getMoment(int j) const {
            if (j == 1) { return getExpectation(); }
            MpfiWrapper den(1.0 + (double) j);

            MpfiWrapper a = range.upper();
            MpfiWrapper b = range.lower();
            den = den * (b-a);
            MpfiWrapper v = pow(b, j+1) - pow(a, j+1);
            MpfiWrapper retVal = 1.0/den * v;
            return retVal;
        }



    };


};

#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_DISTRIBUTIONINFO_HH
