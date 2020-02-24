//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_STATEABSTRACTION_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_STATEABSTRACTION_HH

#include "MultivariatePoly.hh"

namespace PolynomialForms{

    typedef std::shared_ptr< std::map<int, DistributionInfoPtr > > NoiseSymbDistributionInfo;

    class StateAbstraction {
    protected:
        std::map<int, MultivariatePoly> stateVarMap;
        std::map<int, DistributionInfoPtr >  distributionInfo;
        int numNoiseSymbols;
        int maxDegree;

    public:
        StateAbstraction(int maxDegree_=4): numNoiseSymbols(0), maxDegree(maxDegree_) {};

        int getMaxDegree(){
            return maxDegree;
        }

        int getNewNoiseSymbol(){
            int id = numNoiseSymbols;
            numNoiseSymbols ++;
            return id;
        }

        std::map<int, DistributionInfoPtr > const& getNoiseSymbolInfoMap() const {
            return distributionInfo;
        }

        void initialize(std::map<int, DistributionInfoPtr> const & distribMap ){
            for (auto p: distribMap){
                int varId = p.first;
                DistributionInfoPtr dPtr = p.second;
                int noiseSymbID = getNewNoiseSymbol();
                MultivariatePoly mp( MpfiWrapper(1.0), noiseSymbID);
                stateVarMap.insert(make_pair(varId, mp));
                distributionInfo.insert(make_pair(noiseSymbID, dPtr));
            }
            return;
        }

        void addNewNoiseSymbol(int varId, DistributionInfoPtr  dPtr){
            distributionInfo.insert(make_pair(varId, dPtr));
        }

        MultivariatePoly getPolynomialForVar(int varID){
            auto it = stateVarMap.find(varID);
            assert(  it != stateVarMap.end());
            return it -> second;
        }

        std::map<int, MpfiWrapper> getRangeMapForNoiseSymbols(){
            std::map<int, MpfiWrapper> retMap;
            for (auto p: distributionInfo){
                retMap.insert(make_pair(p.first, p.second->getRange()));
            }
            return retMap;
        }

        void setStateMap(std::map<int, MultivariatePoly> const & newStateMap){
            stateVarMap.clear();
            stateVarMap = newStateMap;
        }
    };

    typedef std::shared_ptr<StateAbstraction> StateAbstractionPtr;

};

#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_STATEABSTRACTION_HH
