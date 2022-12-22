#ifndef CONSTANTS
#define CONSTANTS


#include <complex>


// The extension to be added to .dat files to indicate which
// order of accuracy is used for the simulation.
const std::string EXTENSION = "";
//const std::string EXTENSION = "_4";
//const std::string EXTENSION = "_5";


// Constants are defined using SI units.
const double GAMMA_H = 2.675222e8; // gyromagnetic ratio of 1H in rad/(s T)
const double H_BAR = 1.054571e-34; // hbar
const double TWO_PI = 2 * 3.141592653; // 2pi
const double MU_0 = 2 * TWO_PI * 1.0e-7; // permeability of free space
const double PI = 3.141592653;         // pi


// Complex numbers that are often used
const std::complex<double> complex0(0.0, 0.0);
const std::complex<double> complex1(1.0, 0.0);
const std::complex<double> complexi(0.0, 1.0);

#endif
