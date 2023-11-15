/*
 * constants.h
 *
 *  Created on: Feb 9, 2023
 *      Author: pemueller
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <math.h>
#include <complex>

const double symmetryPrecision = 5e-4;
const complex<double> complexZero(0,0);
const complex<double> complexI(0,1);
const complex<double> TWOPITimesI = 2 * M_PI * complexI;
const complex<double> PiTimesI = M_PI * complexI;
const double TWOPI = 2 * M_PI;
const int infinityRotationSteps = 4; // max angular momentum of orbitals; minimum 4, because D2h point group would not be fun



#endif /* CONSTANTS_H_ */
