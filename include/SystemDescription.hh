//
// Created by Sriram Sankaranarayanan on 2/13/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_SYSTEMDESCRIPTION_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_SYSTEMDESCRIPTION_HH

#include "MultivariatePoly.hh"
#include "ExprAST.hh"
#include "DistributionInfo.hh"
#include "Queries.hh"
#include "StateAbstraction.hh"

namespace PolynomialForms {

  class StochasticSystem {
  protected:
      int numStateVars;
      std::map<int, string> varNames;
      std::map<string, int> varIDs;
      std::map<int, ExprPtr> updates;
      std::map<int, DistributionInfoPtr> initialDistrib;
      std::vector<Query> queries;

  public:

      StochasticSystem(){}

      /*StochasticSystem(int nVars,
                       std::map<int, string> const & varNames_,
                       std::map<int, ExprPtr> const & updates_,
                       std::map<int, DistributionInfoPtr> const & initialDistrib_):numStateVars(nVars),
                                                                                   varNames(varNames_),
                                                                                   updates(updates_),
                                                                                   initialDistrib(initialDistrib_)
                                                                                   {};*/

      int addVar(std::string varName, DistributionInfoPtr initialDistrib);
      void addUpdate(int varID, ExprPtr rhs);
      int getVarIDFromName(std::string varName);
      void addQuery(Query q){
          queries.push_back(q);
      }

      std::map<int, DistributionInfoPtr > const & getInitialMap() const {
          return initialDistrib;
      }

      StateAbstractionPtr initialize(int maxDegree =4);
      void computeOneStep(StateAbstractionPtr st, int maxDegree = -1);

      void prettyPrintStateAbstraction(std::ostream & what, StateAbstractionPtr st);

      void evaluateQueries(StateAbstractionPtr st);
  };

  typedef std::shared_ptr<StochasticSystem> StochasticSystemPtr;
};
#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_SYSTEMDESCRIPTION_HH
