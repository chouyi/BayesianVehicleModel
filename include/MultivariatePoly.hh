//
// Created by Sriram Sankaranarayanan on 2/3/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_MULTIVARIATEPOLY_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_MULTIVARIATEPOLY_HH

#include <vector>
#include <map>
#include <set>
#include "MpfiWrapper.hh"
#include <iostream>
#include "DistributionInfo.hh"
#include <sstream>

namespace PolynomialForms{

    class PowerProduct {
    protected:
        std::map<int, int> pp;
        std::set<int> indices;

        std::set<int> combineVars(PowerProduct const & other ) const {
            std::set<int> allIndices;
            std::set_union(indices.begin(), indices.end(),
                           other.indices.begin(), other.indices.end(),
                           std::inserter(allIndices, allIndices.begin()));
            return allIndices;
        }

    public:
        PowerProduct() {};
        PowerProduct( PowerProduct const & q):  pp(q.pp), indices(q.indices){}
        int operator() (int id) const {
            auto it = pp.find(id);
            if (it != pp.end())
                return it -> second;
            else
                return 0;
        }

        set<int> const & getIndices() const { return indices; }


        bool isZeroPP() const {
            return totalDegree() == 0;
        }

        /* Check if any of the vars in the power product are in the set vars */
        bool hasIntersectionWith(std::set<int> const & vars ) const {
            return std::any_of(indices.begin(), indices.end(), [&](int i) {
               return vars.find(i) != vars.end();
            } );
        }



        int totalDegree() const {
            int sum = 0;
            for (auto p: pp) {
                sum += p.second;
            }
            return sum;
        }

        bool operator< (PowerProduct const & other) const {
            std::set<int> allIndices = combineVars(other);

            for (int j : allIndices){
                if ((*this)(j) > other(j)){
                    return false;
                }
                if ((*this)(j) < other(j)){
                    return true;
                }
            }

            return false;
        }

        bool operator== (PowerProduct const & other) const {
            std::set<int> allIndices = combineVars(other);
            for (int j : allIndices){
                if ((*this)(j) != other(j)){
                    return false;
                }
            }

            return true;
        }

        void setPower(int varID, int pow){
            assert(pow >= 0);
            if (pow == 0) {
                std::cout << "Warning: attempting to set to zero power" << std::endl;
                return;
            }
            pp[varID] = pow;
            indices.insert(varID);
        }


        PowerProduct multiply(PowerProduct const & other) const {
            PowerProduct retP;
            std::set<int> allIndices = combineVars(other);
            for (int j : allIndices) {
                retP.setPower(j, (*this)(j) + other(j));
            }
            return retP;
        }

        MpfiWrapper evaluate(std::map<int, MpfiWrapper> const & env) const {
            MpfiWrapper retValue(1.0);
            for (auto const p: pp){
                auto it = env.find(p.first);
                assert(it != env.end());
                retValue = retValue * pow(it -> second, p.second);
            }
            return retValue;
        }

        MpfiWrapper expectation(std::map<int, DistributionInfoPtr> const & env) const {
            MpfiWrapper retValue(1.0);
            for (auto const p: pp){
                auto it = env.find(p.first);
                assert(it != env.end());
                retValue = retValue * (it -> second) -> getMoment(p.second);
            }
            return retValue;
        }

        void prettyPrint(ostream & out, std::map<int, string> const & name_env) const {
            std::string sep = "";
            for (auto const p: pp){
                auto it = name_env.find(p.first);
                stringstream ss;
                string name;
                if (it == name_env.end()){
                    ss << "w"<<p.first ;
                    name = ss.str();
                }  else {
                    name = it -> second;
                }

                assert(p.second > 0);
                out << sep;
                out << name ;
                if (p.second > 1){
                    out << "^" <<p.second;
                }
                sep = "* ";
            }
        }





    };
    typedef enum {PLUS_SINE, PLUS_COSINE, MINUS_SINE, MINUS_COSINE} trig_deriv_t;
    class MultivariatePoly {
    protected:
        MpfiWrapper constIntvl;
        std::map<PowerProduct, MpfiWrapper> terms;

        MultivariatePoly computeTrig(
                std::map<int, MpfiWrapper> const &var_env,
                trig_deriv_t derivState
        ) const;
    public:

        MultivariatePoly(){};
        MultivariatePoly(MpfiWrapper const & what):constIntvl(what){};

        MultivariatePoly(MpfiWrapper const & what, int varID){
            PowerProduct pp;
            pp.setPower(varID, 1);
            terms.insert(std::make_pair(pp, what));
        }

        int degree() const {
            int maxDegree = 0;
            for (auto p: terms){
                int d = p.first.totalDegree();
                if (d > maxDegree)
                    maxDegree = d;
            }
            return maxDegree;
        }

        std::map<PowerProduct, MpfiWrapper> getTerms() const { return terms;}
        MpfiWrapper getConstIntvl() const {return constIntvl; }
        std::set<PowerProduct> getPowerProducts() const {
            set<PowerProduct> retSet;
            for (auto p: terms){
                retSet.insert(p.first);
            }
            return retSet;
        }

        void setConst(MpfiWrapper const & coeff){
            constIntvl = coeff;
        }

        void setTerm(PowerProduct const & pp, MpfiWrapper const & coeff){
            if (pp.isZeroPP()){
                constIntvl = coeff;
            } else {
                terms[pp] = coeff;
            }
        }

        void addToConst(MpfiWrapper const & coeff){
            constIntvl = constIntvl + coeff;
        }

        void addToTerm(PowerProduct const & pp, MpfiWrapper const & coeff){
            if (pp.isZeroPP()){
                std::cout << "Warning: addToTerm - zero monomial!" << std::endl;
                constIntvl = constIntvl + coeff;
            } else {
                auto it = terms.find(pp);
                if (it != terms.end()){
                    setTerm(pp, it -> second + coeff);
                } else {
                    setTerm(pp, coeff);
                }
            }
        }

        void scaleAndAddAssign(MpfiWrapper const & scale, MultivariatePoly const & mp, bool handleConst = true){
            for (auto p: mp.terms){
                addToTerm(p.first, p.second * scale);
            }
            // Let us take care of the constant term.
            if (handleConst) {
                constIntvl = constIntvl + scale * mp.constIntvl;
            }
        }

        void scaleAssign(MpfiWrapper const & scale){
            if (scale.upper() <= 0.0 && scale.lower() >= 0.0){
                constIntvl = 0.0;
                terms.clear();
                return;
            } else {
                constIntvl = scale * constIntvl;
                for (auto & p:terms){
                    p.second.set(p.second * scale);
                }
            }
        }

        void multiplyAndAddAssign(MpfiWrapper const & c, PowerProduct const & pp,
                                MultivariatePoly const & mp, bool handleConst = true) {
            if (pp.isZeroPP()){
                std::cout << " Multiply and add assign with constant! " << std::endl;
                scaleAndAddAssign(c, mp, handleConst);
            } else {
                for (auto p: mp.terms) {
                    PowerProduct prod = pp.multiply(p.first);
                    addToTerm(prod, p.second * c);
                }
                // Take care of the constant term
                if (handleConst) {
                    addToTerm(pp, constIntvl * c);
                }
            }
        }
        // Evaluate the polynomial under an environment
        MpfiWrapper evaluate(std::map<int, MpfiWrapper> const & var_env) const;
        MpfiWrapper expectation(std::map<int, DistributionInfoPtr> const & env) const;
        MultivariatePoly multiply (MultivariatePoly const & what) const;
        MultivariatePoly rangeMultiply(MultivariatePoly const & what, std::map<int, MpfiWrapper> const & var_env) const;
        MultivariatePoly squarePoly() const ;
        MultivariatePoly powPoly(int j) const;
        MultivariatePoly sine(std::map<int, MpfiWrapper> const &var_env) const ;
        MultivariatePoly cosine(std::map<int, MpfiWrapper> const &var_env) const;
      //  MultivariatePoly exp(std::map<int, MpfiWrapper> const &var_env) const ;
        void truncateAssign(int maxDegree,  std::map<int, MpfiWrapper> const &var_env );
        MultivariatePoly truncate(int maxDegree, std::map<int, MpfiWrapper> const &var_env) const;
        void centerAssign(std::map<int, MpfiWrapper> const &var_env);
        void prettyPrint(ostream & out, std::map<int, string> const & name_env) const;
        MpfiWrapper constantIntvl() const { return constIntvl; };


        MultivariatePoly reciprocal(std::map<int, MpfiWrapper> const & var_env) const;

    };
};

#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_MULTIVARIATEPOLY_HH
