#ifndef TENSOR_CLASS
#define TENSOR_CLASS


#include <cmath>


// This files contains a R2_IST class, which describes a rank 2 irreducible 
// spherical tensor and has member functions to perform rotations.
// A Euler class and Wigner class are also defined, to store euler angles and 
// to calculate the reduced wigner d-functions respectively.


// Set of euler (angles alpha, beta, gamma) that can be used to desribe 
// rotations. The conventions of Rose are used.
class Euler {

public:

    // Euler angles
    double alpha;
    double beta;
    double gamma;

    Euler();            // Default constructor
    Euler(double a, double b, double c); // Constructor
    Euler(const Euler &e);    // Copy constructor
    ~Euler();            // Destructor


};

// Reduced wigner functions for the rotation of irreducible spherical tensors
// of rank 2. The constructor is given a value of beta and the d_mm' are 
// calculated and stored in variables of type double
class Wigner {

public:

    // Elements of the wigner reduced functions
    double d_22, d_12, d_02, d_m12, d_m22;
    double d_21, d_11, d_01, d_m11, d_m21;
    double d_20, d_10, d_00, d_m10, d_m20;

    // The d_xm1 and d_xm2 are not necessary, as the A_2m1 and A_2m2
    // components of R2_IST can be obtained from symmetry relations.

    Wigner();            // Default constructor
    Wigner(double beta);        // Constructor
    Wigner(const Wigner &w);    // Copy constructor
    ~Wigner();            // Destructor

};


class R2_IST {

public:

    // The components are stored in the followong order:
    // A_22, A_21, A_20, A_2m1, A_2m2
    std::complex<double> components[5];


    R2_IST();            // default constructor
    R2_IST(const R2_IST &tensor, const Euler &euler); // Constructor that takes
    // as an input
    // an initital R2_IST and a set of Euler angles, and that creates a new
    // R2_IST by rotation of the initial one
    R2_IST(const R2_IST &tensor);    // Copy constructor
    ~R2_IST();    // Destructor

};

#endif
