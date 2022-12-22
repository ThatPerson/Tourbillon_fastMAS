#include "constants.h"
#include <cmath>
#include <complex>


using namespace std;


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


// Default constructor
Euler::Euler() {

    alpha = 0.0;
    beta = 0.0;
    gamma = 0.0;

}


// Constructor
Euler::Euler(double a, double b, double c) {

    alpha = a;
    beta = b;
    gamma = c;

}


// Copy constructor
Euler::Euler(const Euler &e) {

    alpha = e.alpha;
    beta = e.beta;
    gamma = e.gamma;

}


// Destructor
Euler::~Euler() {

}


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


// Default constructor
Wigner::Wigner() {

    d_22 = 0.0;
    d_12 = 0.0;
    d_02 = 0.0;
    d_m12 = 0.0;
    d_m22 = 0.0;

    d_21 = 0.0;
    d_11 = 0.0;
    d_01 = 0.0;
    d_m11 = 0.0;
    d_m21 = 0.0;

    d_20 = 0.0;
    d_10 = 0.0;
    d_00 = 0.0;
    d_m10 = 0.0;
    d_m20 = 0.0;

}


// Constructor
Wigner::Wigner(double beta) {

    // calculates the relevant trigonometric functions
    double cos_b = cos(beta);
    double sin_b = sin(beta);
    double sq_cos_b = cos_b * cos_b;
    double sq_sin_b = sin_b * sin_b;
    double sin_2b = sin(2 * beta);


    // Only nine elements need to be calculated explicitely, the others can
    // be obtained from symmetry relations
    d_22 = 0.25 * (1 + cos_b) * (1 + cos_b);
    d_12 = 0.5 * (1 + cos_b) * sin_b;
    d_02 = sqrt(3.0 / 8.0) * sq_sin_b;
    d_m12 = 0.5 * (1 - cos_b) * sin_b;
    d_m22 = 0.25 * (1 - cos_b) * (1 - cos_b);

    d_21 = -1.0 * d_12;
    d_11 = 0.5 * (cos_b - 1) + sq_cos_b;
    d_01 = sqrt(3.0 / 8.0) * sin_2b;
    d_m11 = 0.5 * (1 + cos_b) - sq_cos_b;
    d_m21 = d_m12;

    d_20 = d_02;
    d_10 = -1.0 * d_01;
    d_00 = 0.5 * (3 * sq_cos_b - 1);
    d_m10 = d_01;
    d_m20 = d_02;

}


// Copy constructor
Wigner::Wigner(const Wigner &w) {

    d_22 = w.d_22;
    d_12 = w.d_12;
    d_02 = w.d_02;
    d_m12 = w.d_m12;
    d_m22 = w.d_m22;

    d_21 = w.d_21;
    d_11 = w.d_11;
    d_01 = w.d_01;
    d_m11 = w.d_m11;
    d_m21 = w.d_m21;

    d_20 = w.d_20;
    d_10 = w.d_10;
    d_00 = w.d_00;
    d_m10 = w.d_m10;
    d_m20 = w.d_m20;

}


// Destructor
Wigner::~Wigner() {

}


class R2_IST {

public:

    // The components are stored in the following order:
    // A_22, A_21, A_20, A_2m1, A_2m2
    std::complex<double> components[5];


    R2_IST();            // default constructor
    R2_IST(const R2_IST &tensor, const Euler &euler); // Constructor that takes
    // as an input an initital R2_IST and a set of Euler angles, and that
    // creates a new R2_IST by rotation of the initial one
    R2_IST(const R2_IST &tensor);    // Copy constructor
    ~R2_IST();    // Destructor

};


// Default constructor
R2_IST::R2_IST() {

    //  components = new complex<double> [5];

    for (int i = 0; i < 5; i++)
        components[i] = complex0;

}


// copy construcor
R2_IST::R2_IST(const R2_IST &tensor) {

    //  components = new complex<double> [5];

    for (int i = 0; i < 5; i++)
        components[i] = tensor.components[i];

}


// Constructor that performs a rotation
R2_IST::R2_IST(const R2_IST &tensor, const Euler &euler) {

    //  components = new complex<double> [5];

    Wigner wigner(euler.beta);

    std::complex<double> A_22 = tensor.components[0]
                                * exp(complexi * -2.0 * euler.alpha);
    std::complex<double> A_21 = tensor.components[1]
                                * exp(complexi * -1.0 * euler.alpha);
    std::complex<double> A_20 = tensor.components[2];
    // Symmetry properties of irreducible tensors can be used
    std::complex<double> A_2m1 = -1.0 * conj(A_21);
    std::complex<double> A_2m2 = conj(A_22);

    components[0] = exp(complexi * -2.0 * euler.gamma) *
                    (A_22 * wigner.d_22 + A_21 * wigner.d_12 + A_20 * wigner.d_02 +
                     A_2m1 * wigner.d_m12 + A_2m2 * wigner.d_m22);

    components[1] = exp(complexi * -1.0 * euler.gamma) *
                    (A_22 * wigner.d_21 + A_21 * wigner.d_11 + A_20 * wigner.d_01 +
                     A_2m1 * wigner.d_m11 + A_2m2 * wigner.d_m21);

    components[2] =
            (A_22 * wigner.d_20 + A_21 * wigner.d_10 + A_20 * wigner.d_00 +
             A_2m1 * wigner.d_m10 + A_2m2 * wigner.d_m20);

    // Symmetry properties of irreducible tensors can be used
    components[3] = -1.0 * conj(components[1]);
    components[4] = conj(components[0]);

}


R2_IST::~R2_IST() {

    //  delete[] components;

}
