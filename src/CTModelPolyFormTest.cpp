//
// Created by Sriram Sankaranarayanan on 2/24/20.
//

#include "CTModelPolyForm.hh"

int main(){
    CTModelPolynomialForm ctmp(0.2, 12.1, 12.3, -0.1, 0.1, 0.4, 10);
    ctmp.computeNStepForm(true);
}