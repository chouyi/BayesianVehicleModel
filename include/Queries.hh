//
// Created by Sriram Sankaranarayanan on 2/14/20.
//

#ifndef POLYNOMIALFORMUNCERTAINTYPROPAGATION_QUERIES_HH
#define POLYNOMIALFORMUNCERTAINTYPROPAGATION_QUERIES_HH

#include <string>
#include "ExprAST.hh"

namespace PolynomialForms {
    typedef enum {PROB_QUERY, EXPECT_QUERY } query_type_t;
    class Query {
    protected:
        query_type_t qType;
        std::string id;
        ExprPtr queryExpr; // we seek its expectation

    public:
        Query(query_type_t qType_, std::string qid, ExprPtr e_) : qType(qType_), id(qid), queryExpr(e_) {}
        ExprPtr getExpr(){ return queryExpr; }
        string getID() const { return id;}
        query_type_t  getType() const { return qType; }
    };

    inline Query probabilityQuery(std::string qid, ExprPtr e_) {
        return Query(PROB_QUERY, qid, e_);
    }


    inline Query expectationQuery(std::string qid, ExprPtr e_) {
        return Query(EXPECT_QUERY, qid, e_);
    }


};
#endif //POLYNOMIALFORMUNCERTAINTYPROPAGATION_QUERIES_HH
