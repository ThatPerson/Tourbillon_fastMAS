// The following preprocessor instructions can be used to include in the
// reduced Liouville spaceproduct-operators that involve five single-spin 
// operators. 

#include "constants.h"
#include "hamiltonian.h"
#include "crystal.h"
#include <iostream>
#include <cmath>
#include <complex>
#include <algorithm>
#include <vector>
#include <bitset>
#include <iomanip>
#include <fstream>
#include <unordered_map>

using namespace std;

struct permijk {
    int *i;
    int *j;
    int *k;
};
struct permijkl {
    int *i;
    int *j;
    int *k;
    int *l;
};
// A Density_matrix stores, in a series of multidimensional arrays, the 
// coefficient b_r of the expansion of the density matrix in a basis set of
// product operators: sigma = sum_r( b_r * B_r ). 
// The following notation is used for the product operators: 
// z : I_z, p : I_+ * (1/sqrt(2)), m = I_- * (1/sqrt(2))
// The product operators have a 2^(q-1) normalisation factor, where q is the 
// number of spins; 
// e.g. zpm[i][j][k] stand for 4 * I_iz * (I_j+ / sqrt(2)) * (I_k- / sqrt(2))

// A Density_matrix has an evolve member function, that perfoms a numerical 
// integration of the LVN equation using a Suzuki-Trotter scheme
class Density_matrix {

public:

    /* Here, the original Tourbillon matrices are commented out. As in Perras 2019,
       an unordered map has been implemented. The index in the map is derived 
       from a list of the interacting spins using the private 'get_hash()' function*/
    //std::complex<double> *z;    // Element [i] is the coefficient of I_iz
    std::unordered_map<int, std::complex<double>> z;

    //std::complex<double> **pm;       // Element [i][j] is the coefficient of
    std::unordered_map<int, std::complex<double>> pm;
    // I_i+ * I_j-

    //std::complex<double> **zz;    // Element [i][j] is the coefficient of
    std::unordered_map<int, std::complex<double>> zz;
    // 2 * I_iz * I_jz

    std::unordered_map<int, std::complex<double>> zpm;    // Element [i][j][k] is the coefficient of
    // 2 * I_iz * I_j+ * I_k-

    std::unordered_map<int, std::complex<double>> zzz;    // Element [i][j][k] is the coefficient of
    // 4 * I_iz * I_jz * I_kz.

    std::unordered_map<int, std::complex<double>> ppmm;    // Element [i][j][k][l] is the
    //coefficient of 2 * I_i+ * I_j+ * I_k- * I_l-

    std::unordered_map<int, std::complex<double>> zzpm;    // Element [i][j][k][l] is the
    // coefficient of 4 * I_iz * I_jz * I_k+ * I_l-
    std::vector<std::complex<double> *> saturate_z;

    //std::complex<double> ****zzzz;    // Element [i][j][k][l] is the
    // coefficient of 8 * I_iz * I_jz * I_kz * I_lz



    // Constructor function
    Density_matrix(int number_of_nuclei, Crystal &crystal);
	/* Should take Crystal &crystal */

    // Performs the numerical integration of the LVN equation with a
    // Suzuki-Trotter algorithm of order 1 or 2.
    void evolve_ST1(const Hamiltonian &hamiltonian, double timestep, const Crystal &crystal);
    void evolve_CS(const Crystal &crystal, const Hamiltonian &hamiltonian, double timestep);
    void evolve_ST2(const Hamiltonian &hamiltonian, double timestep, const Crystal &crystal);

    // sets all the elements to zero
    void reset(const Crystal &crystal);
    int number_of_spins;        // Number of spins
    void write_continuation(string filename);
    void read_continuation(string filename);
    void saturate(std::vector<int>& polarised_spins);
    void store_br(double **br_zz, double **br_pm, double **br_mp, const Crystal &crystal);
private:

    int get_hash(std::vector<int> spins);
    // The following functions are called by the main evolve function. update_R
    // calculates the matrix elements for the rotation matrices, and evolve
    // implements their action, for the dipolar Hamiltonian between spins i and j
    void update_RHH(const Hamiltonian &hamiltonian,
                    double timestep,
                    int i,
                    int j);

    void evolve_HH(int i, int j, const Crystal &crystal, const Hamiltonian &hamiltonian, double timestep);


    // The evolution makes repeated use of the same matrix elements
    double RHHzqc_11, RHHzqc_21;
    std::complex<double> RHHzqc_31, RHHzqc_41;
    double RHHpoqc_11, RHHpoqc_21;
    std::complex<double> RHHpoqc_31, RHHpoqc_41;
    double RHHmoqc_11, RHHmoqc_21;
    std::complex<double> RHHmoqc_31, RHHmoqc_41;

    // The following functions multiply the vector X = [a;b;c;d]
    // by a 4 by 4 rotation matrix
    void RHHzqcX(std::complex<double> &a, std::complex<double> &b,
                 std::complex<double> &c, std::complex<double> &d);

    void RHHpoqcX(std::complex<double> &a, std::complex<double> &b,
                  std::complex<double> &c, std::complex<double> &d);

    void RHHmoqcX(std::complex<double> &a, std::complex<double> &b,
                  std::complex<double> &c, std::complex<double> &d);


};




// Constructor function for the density matrix. It allocates the memory
// for the arrays containing the coefficients of the expansion of the 
// density matrix. 
// Only the number of spins needs to be known to contruct a density matrix.

/* : Needs to know crystal dimensions */
Density_matrix::Density_matrix(int number_of_nuclei, Crystal &crystal) {

    number_of_spins = number_of_nuclei;

    // Allocate memory for the arrays storing the coefficients for the expansion
    // of the density matrix. Each array has a shape that avoids redundancies.



    // Initialise all the elements to zero
    reset(crystal);
    // Initialise the rotation matrices to zero
    RHHzqc_11 = 0.0;
    RHHzqc_21 = 0.0;
    RHHzqc_31 = complex0;
    RHHzqc_41 = complex0;
    RHHpoqc_11 = 0.0;
    RHHpoqc_21 = 0.0;
    RHHpoqc_31 = complex0;
    RHHpoqc_41 = complex0;
    RHHmoqc_11 = 0.0;
    RHHmoqc_21 = 0.0;
    RHHmoqc_31 = complex0;
    RHHmoqc_41 = complex0;


} // constructor

inline int Density_matrix::get_hash(std::vector<int> spins) {
    int c = 1;
    int hash;
    for (auto spin : spins) {
        hash += c * spin;
        c *= number_of_spins;
    }
    return hash;
}

void Density_matrix::write_continuation(string filename) {
    // most definitely not the best way to do it.
    std::vector<pair<string, unordered_map<int, complex<double>>*>> q;
    q.emplace_back("z", &z);
    q.emplace_back("pm", &pm);
    q.emplace_back("zz", &zz);
    q.emplace_back("zpm", &zpm);
    q.emplace_back("zzz", &zzz);
    q.emplace_back("ppmm", &ppmm);
    q.emplace_back("zzpm", &zzpm);

    for (auto qp : q) {
        FILE *f = fopen((filename + qp.first + ".cont").c_str(), "wb");
        for (auto key: *(qp.second)) {
            fwrite(&(key.first), sizeof(int), 1, f);
            fwrite(&(key.second), sizeof(double), 1, f);
        }
        fclose(f);
    }

   /* FILE *f = fopen("map", "wb");
    for(iter=map.begin(); iter!=map.end(); iter++){
        fwrite(&(iter->first), 8, 1, f);
        fwrite(&(iter->second), 1, 1, f);
    }
    fclose(f);*/
}

void Density_matrix::read_continuation(string filename) {
    /* Reading and writing continuations does not work well as it doesn't preserve the state of the 
       propagator - so the evolution post reading a continuation will not be the same as it was before.
       Do not use!
     */
    std::vector<pair<string, unordered_map<int, complex<double>>*>> q;
    q.emplace_back("z", &z);
    q.emplace_back("pm", &pm);
    q.emplace_back("zz", &zz);
    q.emplace_back("zpm", &zpm);
    q.emplace_back("zzz", &zzz);
    q.emplace_back("ppmm", &ppmm);
    q.emplace_back("zzpm", &zzpm);

    for (auto qp : q) {
        FILE *f = fopen((filename + qp.first + ".cont").c_str(), "rb");
        int key;
        double val;
        while (fread(&key, sizeof(int), 1, f)) {
            fread(&val, sizeof(double), 1, f);
            //*(qp.second)[key] = val;
            if (qp.first == "z") { z[key] = val; }
            else if (qp.first == "pm") { pm[key] = val; }
            else if (qp.first == "zz") { zz[key] = val; }
            else if (qp.first == "zpm") { zpm[key] = val; }
            else if (qp.first == "zzz") { zzz[key] = val; }
            else if (qp.first == "ppmm") { ppmm[key] = val; }
            else if (qp.first == "zzpm") { zzpm[key] = val; }
        }
        fclose(f);
    }
}

void Density_matrix::reset(const Crystal &crystal) {


    // sets all the elements to zero


    for (int i = 0; i < number_of_spins; i++) {
        z[i] = complex0;
        for (auto j : crystal.nearest_neighbours[i]) {
            pm[get_hash({i,j})] = complex0;
            pm[get_hash({j,i})] = complex0;
            zz[get_hash({j, i})] = complex0;
        }
    }

    for (int i = 0; i < number_of_spins; i++) {
        for (auto j : crystal.nearest_neighbours[i]) {
            std::set<int> kset;
            std::set_intersection(crystal.nearest_neighbours[i].begin(), crystal.nearest_neighbours[i].end(),
                                  crystal.nearest_neighbours[j].begin(), crystal.nearest_neighbours[j].end(),
                                  std::inserter(kset, kset.begin()));
            for (auto k : kset) {
                if ((k <= j) || (k <= i))
                    continue;

                zpm[get_hash({i,j,k})] = complex0;
                zpm[get_hash({i,k,j})] = complex0;
                zpm[get_hash({j,i,k})] = complex0;
                zpm[get_hash({j,k,i})] = complex0;
                zpm[get_hash({k,i,j})] = complex0;
                zpm[get_hash({k,j,i})] = complex0;
                zzz[get_hash({k,j,i})] = complex0; // k > j > i
            }
        }
    }

    std::vector<permijkl> permutations;
    int a, b, c, d;
    permutations.push_back({&a, &b, &c, &d});
    permutations.push_back({&a, &b, &d, &c});
    permutations.push_back({&a, &c, &b, &d});
    permutations.push_back({&a, &c, &d, &b});
    permutations.push_back({&a, &d, &b, &c});
    permutations.push_back({&a, &d, &c, &b});
    permutations.push_back({&c, &a, &b, &d});
    permutations.push_back({&c, &a, &d, &b});
    permutations.push_back({&c, &b, &a, &d});
    permutations.push_back({&c, &b, &d, &a});
    permutations.push_back({&c, &d, &a, &b});
    permutations.push_back({&c, &d, &b, &a});
    permutations.push_back({&b, &a, &c, &d});
    permutations.push_back({&b, &a, &d, &c});
    permutations.push_back({&b, &d, &a, &c});
    permutations.push_back({&b, &d, &c, &a});
    permutations.push_back({&b, &c, &a, &d});
    permutations.push_back({&b, &c, &d, &a});
    permutations.push_back({&d, &a, &c, &b});
    permutations.push_back({&d, &a, &b, &c});
    permutations.push_back({&d, &b, &a, &c});
    permutations.push_back({&d, &b, &c, &a});
    permutations.push_back({&d, &c, &a, &b});
    permutations.push_back({&d, &c, &b, &a});

    for (int i = 0; i < number_of_spins; i++) {
        for (auto j : crystal.nearest_neighbours[i]) {
            std::set<int> kset;
            std::set_intersection(crystal.nearest_neighbours[i].begin(), crystal.nearest_neighbours[i].end(),
                                  crystal.nearest_neighbours[j].begin(), crystal.nearest_neighbours[j].end(),
                                  std::inserter(kset, kset.begin()));
            for (auto k : kset) {
                if ((k <= i) || (k <= j)) {
                    continue;
                }
                std::set<int> lset;
                std::set_intersection(kset.begin(), kset.end(),
                                      crystal.nearest_neighbours[k].begin(), crystal.nearest_neighbours[k].end(),
                                      std::inserter(lset, lset.begin()));
                for (auto l : lset) {
                    if (l <= k) continue;

                    a = i; b = j; c = k; d = l;
                    for (auto perm : permutations) {
                        if (*(perm.i) > *(perm.j) && *(perm.k) > *(perm.l)) ppmm[get_hash({*(perm.i),*(perm.j),*(perm.k),*(perm.l)})] = complex0;
                        if (*(perm.i) > *(perm.j)) zzpm[get_hash({*(perm.i),*(perm.j),*(perm.k),*(perm.l)})] = complex0;

                    }
                }
            }
        }
    }



}


void Density_matrix::evolve_CS(const Crystal &crystal, const Hamiltonian &hamiltonian, double timestep) {
    auto *compexp_p = new complex<double>[number_of_spins];
    auto *compexp_m = new complex<double>[number_of_spins];
    double shifts_Hz;
    for (int i = 0 ; i < number_of_spins; i++) {
        //hamiltonian.omega_D_L_20[i][j]
        double cs_total = crystal.csiso_ppm[i] + hamiltonian.omega_CSA_L_20[i]; // calculate instantaneous chemical shift
        
        //cs_total += hamiltonian.omega_CSA_L_2m2[i] + hamiltonian.omega_CSA_L_2p2[i];

        shifts_Hz = cs_total * crystal.field_strength;
       // printf("Shifts %f\n", shifts_Hz);
        compexp_p[i] = exp(complexi * 2. * M_PI * shifts_Hz * timestep / 1000.);
        compexp_m[i] = exp(-complexi * 2. * M_PI * shifts_Hz * timestep / 1000.);
    }

    std::vector<permijkl> permutations;
    permutations.reserve(24);
    int a, b, c, d;
    permutations.push_back({&a, &b, &c, &d});
    permutations.push_back({&a, &b, &d, &c});
    permutations.push_back({&a, &c, &b, &d});
    permutations.push_back({&a, &c, &d, &b});
    permutations.push_back({&a, &d, &b, &c});
    permutations.push_back({&a, &d, &c, &b});
    permutations.push_back({&c, &a, &b, &d});
    permutations.push_back({&c, &a, &d, &b});
    permutations.push_back({&c, &b, &a, &d});
    permutations.push_back({&c, &b, &d, &a});
    permutations.push_back({&c, &d, &a, &b});
    permutations.push_back({&c, &d, &b, &a});
    permutations.push_back({&b, &a, &c, &d});
    permutations.push_back({&b, &a, &d, &c});
    permutations.push_back({&b, &d, &a, &c});
    permutations.push_back({&b, &d, &c, &a});
    permutations.push_back({&b, &c, &a, &d});
    permutations.push_back({&b, &c, &d, &a});
    permutations.push_back({&d, &a, &c, &b});
    permutations.push_back({&d, &a, &b, &c});
    permutations.push_back({&d, &b, &a, &c});
    permutations.push_back({&d, &b, &c, &a});
    permutations.push_back({&d, &c, &a, &b});
    permutations.push_back({&d, &c, &b, &a});

    /* See reset() - need to do the same here! */
    // For pm[i][j]
    for (int i = 0; i < number_of_spins; i++) {
        z[i] += crystal.Rz * timestep * (complex1 - z[i]);
        for (auto j : crystal.nearest_neighbours[i]) {
            //printf("%d, %d\n", i, j);
            int hash;
            hash = get_hash({i, j});
            pm[hash] *= compexp_p[i] * compexp_m[j];
            pm[hash] += -pm[hash] * timestep * complex1 * 2. * crystal.Rpm;

            hash = get_hash({j, i});
            pm[hash] *= compexp_p[j] * compexp_m[i];
            pm[hash] += -pm[hash] * timestep * complex1 * 2. * crystal.Rpm;

            hash = get_hash({j, i});
            zz[hash] += -zz[hash] * timestep * complex1 * (crystal.Rz + crystal.Rz);

            //pm[get_hash({i, j})] *= compexp_p[i] * compexp_m[j];
            //pm[get_hash({j, i})] *= compexp_p[j] * compexp_m[i];
//        pm[i] = new complex<double>[number_of_spins];

            // now we only want the union of i and j neighbours
            std::set<int> kset;
            std::set_intersection(crystal.nearest_neighbours[i].begin(), crystal.nearest_neighbours[i].end(),
                                  crystal.nearest_neighbours[j].begin(), crystal.nearest_neighbours[j].end(),
                                  std::inserter(kset, kset.begin()));
            for (auto k : kset) {
                /* j set will contain i, i set will contain j.
                 * So k will contain both - but we don't want to
                 * consider self */
                if ((k == j) || (k == i))
                    continue;
                /* Prevent double counting by checking order. */
                if (!((k > j) && (j > i)))
                    continue;

                // k > j > i
                hash = get_hash({k, j, i});
                zzz[hash] += -zzz[hash] * timestep * complex1 * (crystal.Rz + crystal.Rz + crystal.Rz);

                hash = get_hash({i, j, k});
                zpm[hash] *= compexp_p[j] * compexp_m[k];
                zpm[hash] += -zpm[hash] * timestep * complex1 * (crystal.Rz + crystal.Rpm + crystal.Rpm);

                hash = get_hash({i, k, j});
                zpm[hash] *= compexp_p[k] * compexp_m[j];
                zpm[hash] += -zpm[hash] * timestep * complex1 * (crystal.Rz + crystal.Rpm + crystal.Rpm);

                hash = get_hash({k, i, j});
                zpm[hash] *= compexp_p[i] * compexp_m[j];
                zpm[hash] += -zpm[hash] * timestep * complex1 * (crystal.Rz + crystal.Rpm + crystal.Rpm);

                hash = get_hash({k, j, i});
                zpm[hash] *= compexp_p[j] * compexp_m[i];
                zpm[hash] += -zpm[hash] * timestep * complex1 * (crystal.Rz + crystal.Rpm + crystal.Rpm);

                hash = get_hash({j, k, i});
                zpm[hash] *= compexp_p[k] * compexp_m[i];
                zpm[hash] += -zpm[hash] * timestep * complex1 * (crystal.Rz + crystal.Rpm + crystal.Rpm);

                hash = get_hash({j, i, k});
                zpm[hash] *= compexp_p[i] * compexp_m[k];
                zpm[hash] += -zpm[hash] * timestep * complex1 * (crystal.Rz + crystal.Rpm + crystal.Rpm);


                std::set<int> lset;
                std::set_intersection(crystal.nearest_neighbours[k].begin(), crystal.nearest_neighbours[k].end(),
                                      kset.begin(), kset.end(),
                                      std::inserter(lset, lset.begin()));
                for (auto l : lset) {
                    if (l <= k)
                        continue;
                 
                    a = i; b = j; c = k; d = l;
                    for (auto perm : permutations) {
                        if (*(perm.i) > *(perm.j) && *(perm.k) > *(perm.l)) {
                            hash = get_hash({*(perm.i),*(perm.j),*(perm.k),*(perm.l)});
                            ppmm[hash] *= compexp_p[*(perm.i)] * compexp_p[*(perm.j)] * compexp_m[*(perm.k)] * compexp_m[*(perm.l)];
                            ppmm[hash] += -ppmm[hash] * timestep * complex1 * (crystal.Rpm * 4);
                        }
                        if (*(perm.i) > *(perm.j)) {
                            hash = get_hash({*(perm.i),*(perm.j),*(perm.k),*(perm.l)});
                            zzpm[hash] *= compexp_p[*(perm.k)] * compexp_m[*(perm.l)];
                            zzpm[hash] += -zzpm[hash] * timestep * complex1 * (crystal.Rz + crystal.Rz + crystal.Rpm + crystal.Rpm);
                        }
                    }
                }
            }
        }
    }


    //delete[] compexp_p;
    //delete[] compexp_m;
}

// The evolve function performs the numerical integration of the LVN equation
// For this function to be called, a Hamiltonian must be defined, and 
// memory must be allocated for the elements of the density matrix. Also, a 
// timestep must be defined.
// This global "evolve" function calls a series of evolve subfunctions, one
// for each interaction in the total Hamiltonian. The order in which the
// evolve subfunctions are called make it possible to use higher-order
// integration schemes, according to the Suzuki-Trotter algorithm
void Density_matrix::evolve_ST1(const Hamiltonian &hamiltonian,
                                double timestep, const Crystal &crystal) {


	/* : Only loop over the neighbour set of i */
    for (int i = 0; i < number_of_spins; i++) {
        for (auto j : crystal.nearest_neighbours[i]) {
            update_RHH(hamiltonian, timestep, i, j);
            evolve_HH(i, j, crystal, hamiltonian, timestep);
        }
    }

} // evolve_ST1


void Density_matrix::evolve_ST2(const Hamiltonian &hamiltonian,
                                double timestep, const Crystal &crystal) {

    double halftimestep = timestep / 2.0;

	/* Only loop over the neighbour set of i */
    for (int i = 0; i < number_of_spins; i++) {
        for (auto j : crystal.nearest_neighbours[i]) {

            update_RHH(hamiltonian, halftimestep, i, j);
            evolve_HH(i, j, crystal, hamiltonian, halftimestep);
        }
    }

    for (int i = number_of_spins - 1; i >= 0; i--) {
        for (auto j: crystal.nearest_neighbours[i]) {

            update_RHH(hamiltonian, halftimestep, i, j);
            evolve_HH(i, j, crystal, hamiltonian, halftimestep);
        }
    }

} // evolve_ST



void Density_matrix::evolve_HH(int i, int j, const Crystal &crystal, const Hamiltonian &hamiltonian, double timestep) {
    /*auto *compexp_p = new complex<double>[number_of_spins];
    auto *compexp_m = new complex<double>[number_of_spins];
    double shifts_Hz;
    for (int i = 0 ; i < number_of_spins; i++) {
        shifts_Hz = crystal.shifts_ppm[i] * crystal.field_strength;
        // printf("Shifts %f\n", shifts_Hz);
        compexp_p[i] = exp(complexi * 2. * M_PI * shifts_Hz * timestep / 1000.);
        compexp_m[i] = exp(-complexi * 2. * M_PI * shifts_Hz * timestep / 1000.);
    }*/

    // pm
    RHHzqcX(z[i],
            z[j],
            pm[get_hash({i, j})],
            pm[get_hash({j, i})]); // p(i)m(j)


    /* : Replace this with the same set iteration as in reset().
     * Here it might be a little harder as there's weird ik[2] stuff
     * going on, but I think it should be possible to just replace the loops
     */
    // Loop over coherences that involve more than two spins

    std::set<int> kset;
    std::set_intersection(crystal.nearest_neighbours[i].begin(), crystal.nearest_neighbours[i].end(),
                          crystal.nearest_neighbours[j].begin(), crystal.nearest_neighbours[j].end(),
                          std::inserter(kset, kset.begin()));



    int a, b, c, d;
    std::vector<permijk> permutations;
    permutations.reserve(3);
    permutations.push_back({&a, &b, &c});
    permutations.push_back({&a, &c, &b});
    permutations.push_back({&b, &c, &a});

    std::vector<permijkl> permutations4;
    permutations4.reserve(12);
    permutations4.push_back({&a, &b, &c, &d});
    permutations4.push_back({&a, &b, &d, &c});
    permutations4.push_back({&a, &c, &b, &d});
    permutations4.push_back({&a, &c, &d, &b});
    permutations4.push_back({&a, &d, &b, &c});
    permutations4.push_back({&a, &d, &c, &b});
    permutations4.push_back({&b, &c, &a, &d});
    permutations4.push_back({&b, &c, &d, &a});
    permutations4.push_back({&b, &d, &a, &c});
    permutations4.push_back({&b, &d, &c, &a});
    permutations4.push_back({&c, &d, &a, &b});
    permutations4.push_back({&c, &d, &b, &a});


    for (auto k: kset) {
        // All the spin involved in a coherence must be different
        if ((k == i) || (k == j))
            continue;

        a = i; b = j; c = k;

        for (auto perm : permutations) {
            int ik[2] = {*(perm.i), *(perm.k)};
            sort(ik, ik + 2);
            int jk[2] = {*(perm.j), *(perm.k)};
            sort(jk, jk + 2);

            // zpm
            update_RHH(hamiltonian, timestep, *(perm.i), *(perm.j));
            RHHzqcX(zz[get_hash({ik[1], ik[0]})],
                    zz[get_hash({jk[1], jk[0]})],
                    zpm[get_hash({*(perm.k), *(perm.i), *(perm.j)})],
                    zpm[get_hash({*(perm.k), *(perm.j), *(perm.i)})]);        // z(k)p(i)m(j)
            RHHpoqcX(pm[get_hash({*(perm.i), *(perm.k)})],
                     pm[get_hash({*(perm.j), *(perm.k)})],
                     zpm[get_hash({*(perm.j), *(perm.i), *(perm.k)})],
                     zpm[get_hash({*(perm.i), *(perm.j), *(perm.k)})]); // z(i)p(j)m(k)
            RHHmoqcX(pm[get_hash({*(perm.k), *(perm.i)})],
                     pm[get_hash({*(perm.k), *(perm.j)})],
                     zpm[get_hash({*(perm.j), *(perm.k), *(perm.i)})],
                     zpm[get_hash({*(perm.i), *(perm.k), *(perm.j)})]); // z(i)p(k)m(j)
        }


        // Loop over coherences that involve more than three spins
        std::set<int> lset;
        std::set_intersection(crystal.nearest_neighbours[k].begin(), crystal.nearest_neighbours[k].end(),
                              kset.begin(), kset.end(),
                              std::inserter(lset, lset.begin()));

        //for (int l = 0; l < number_of_spins; l++) {
        for (auto l : lset) {
            // All the spins involved in a coherence must be different
            if ((l == i) || (l == j) || (l == k))
                continue;

            a = i; b = j; c = k; d = l;

            for (auto perm : permutations4) {
                int il[2] = {*(perm.i), *(perm.l)};
                sort(il, il + 2);
                int jl[2] = {*(perm.j), *(perm.l)};
                sort(jl, jl + 2);
                int ikl[3] = {*(perm.i), *(perm.k), *(perm.l)};
                sort(ikl, ikl + 3);
                int jkl[3] = {*(perm.j), *(perm.k), *(perm.l)};
                sort(jkl, jkl + 3);
                int ik[2] = {*(perm.i), *(perm.k)};
                sort(ik, ik+2);
                int jk[2] = {*(perm.j), *(perm.k)};
                sort(jk, jk+2);
                update_RHH(hamiltonian, timestep, *(perm.i), *(perm.j));
                // ppmm
                RHHzqcX(zpm[get_hash({*(perm.i), *(perm.k), *(perm.l)})],
                        zpm[get_hash({*(perm.j), *(perm.k), *(perm.l)})],
                        ppmm[get_hash({ik[1],ik[0],jl[1],jl[0]})],
                        ppmm[get_hash({jk[1],jk[0],il[1],il[0]})]); // pp(ik)mm(jl)

                // zzpm
                if (k > l)
                    RHHzqcX(zzz[get_hash({ikl[2],ikl[1],ikl[0]})],
                            zzz[get_hash({jkl[2],jkl[1],jkl[0]})],
                            zzpm[get_hash({*(perm.k),*(perm.l),*(perm.i),*(perm.j)})],
                            zzpm[get_hash({*(perm.k),*(perm.l),*(perm.j),*(perm.i)})]);    // zz(kl)p(i)m(j)
                RHHpoqcX(zpm[get_hash({*(perm.k),*(perm.i),*(perm.l)})],
                         zpm[get_hash({*(perm.k),*(perm.j),*(perm.l)})],
                         zzpm[get_hash({jk[1],jk[0],*(perm.i),*(perm.l)})],
                         zzpm[get_hash({ik[1],ik[0],*(perm.j),*(perm.l)})]); // zz(ik)p(j)m(l)
                RHHmoqcX(zpm[get_hash({*(perm.k),*(perm.l),*(perm.i)})],
                         zpm[get_hash({*(perm.k),*(perm.l),*(perm.j)})],
                         zzpm[get_hash({jk[1],jk[0],*(perm.l),*(perm.i)})],
                         zzpm[get_hash({ik[1],ik[0],*(perm.l),*(perm.j)})]); // zz(ik)p(l)m(j)
            }




        } // l loop

    } // k loop

} // evolve_HH


// This function updates the matrix elements for the rotation matrix that 
// describes the effect of the dipolar interaction between i and j.
void Density_matrix::update_RHH(const Hamiltonian &hamiltonian,
                                double timestep,
                                int i,
                                int j) {

    double theta = hamiltonian.omega_D_L_20[i][j] * timestep;
    double costheta = cos(theta);
    double sintheta = sin(theta);
    double halftheta = theta * 0.5;
    double coshalftheta = cos(halftheta);
    double sinhalftheta = sin(halftheta);

    RHHzqc_11 = 0.5 * (1 + costheta);
    RHHzqc_21 = 0.5 * (1 - costheta);
    RHHzqc_31 = -0.5 * complexi * sintheta;
    RHHzqc_41 = 0.5 * complexi * sintheta;

    RHHpoqc_11 = costheta * coshalftheta;
    RHHpoqc_21 = -1.0 * sintheta * sinhalftheta;
    RHHpoqc_31 = -1.0 * complexi * sintheta * coshalftheta;
    RHHpoqc_41 = -1.0 * complexi * costheta * sinhalftheta;

    RHHmoqc_11 = RHHpoqc_11;
    RHHmoqc_21 = RHHpoqc_21;
    RHHmoqc_31 = -1.0 * RHHpoqc_31;
    RHHmoqc_41 = -1.0 * RHHpoqc_41;

} // update_RHH


void Density_matrix::RHHzqcX(complex<double> &a, complex<double> &b,
                             complex<double> &c, complex<double> &d) {

    complex<double> a0 = a;
    complex<double> b0 = b;
    complex<double> c0 = c;
    complex<double> d0 = d;

    a = RHHzqc_11 * a0;
    a += RHHzqc_21 * b0;
    a += RHHzqc_31 * c0;
    a += RHHzqc_41 * d0;

    b = RHHzqc_21 * a0;
    b += RHHzqc_11 * b0;
    b += RHHzqc_41 * c0;
    b += RHHzqc_31 * d0;

    c = RHHzqc_31 * a0;
    c += RHHzqc_41 * b0;
    c += RHHzqc_11 * c0;
    c += RHHzqc_21 * d0;

    d = RHHzqc_41 * a0;
    d += RHHzqc_31 * b0;
    d += RHHzqc_21 * c0;
    d += RHHzqc_11 * d0;


} // RHHzqcX


void Density_matrix::RHHpoqcX(complex<double> &a, complex<double> &b,
                              complex<double> &c, complex<double> &d) {

    complex<double> a0 = a;
    complex<double> b0 = b;
    complex<double> c0 = c;
    complex<double> d0 = d;

    a = RHHpoqc_11 * a0;
    a += RHHpoqc_21 * b0;
    a += RHHpoqc_31 * c0;
    a += RHHpoqc_41 * d0;

    b = RHHpoqc_21 * a0;
    b += RHHpoqc_11 * b0;
    b += RHHpoqc_41 * c0;
    b += RHHpoqc_31 * d0;

    c = RHHpoqc_31 * a0;
    c += RHHpoqc_41 * b0;
    c += RHHpoqc_11 * c0;
    c += RHHpoqc_21 * d0;

    d = RHHpoqc_41 * a0;
    d += RHHpoqc_31 * b0;
    d += RHHpoqc_21 * c0;
    d += RHHpoqc_11 * d0;


} // RHHpoqcX


void Density_matrix::RHHmoqcX(complex<double> &a, complex<double> &b,
                              complex<double> &c, complex<double> &d) {

    complex<double> a0 = a;
    complex<double> b0 = b;
    complex<double> c0 = c;
    complex<double> d0 = d;

    a = RHHmoqc_11 * a0;
    a += RHHmoqc_21 * b0;
    a += RHHmoqc_31 * c0;
    a += RHHmoqc_41 * d0;

    b = RHHmoqc_21 * a0;
    b += RHHmoqc_11 * b0;
    b += RHHmoqc_41 * c0;
    b += RHHmoqc_31 * d0;

    c = RHHmoqc_31 * a0;
    c += RHHmoqc_41 * b0;
    c += RHHmoqc_11 * c0;
    c += RHHmoqc_21 * d0;

    d = RHHmoqc_41 * a0;
    d += RHHmoqc_31 * b0;
    d += RHHmoqc_21 * c0;
    d += RHHmoqc_11 * d0;

} // RHHmoqcX


void Density_matrix::saturate(std::vector<int>& polarised_spins) {
    // want to set all z operators to 0.


    if (saturate_z.empty()) {
        for (auto i : polarised_spins) {
            saturate_z.push_back(&(z[i]));
            for (int j = 0; j < number_of_spins; j++) {
                std::vector<int> sp = {j, i};
                sort(sp.begin(), sp.end());
                saturate_z.push_back(&(zz[get_hash(sp)]));
                for (int k = 0; k < number_of_spins; k++) {
                    sp.push_back(k);
                    sort(sp.begin(), sp.end());
                    saturate_z.push_back(&(zzz[get_hash(sp)]));
                    saturate_z.push_back(&(zpm[get_hash({i, j, k})]));
                    saturate_z.push_back(&(zpm[get_hash({i, k, j})]));
                    for (int l = 0; l < number_of_spins; l++) {
                        int ij[2] = {i, j};
                        int ik[2] = {i, k};
                        int il[2] = {i, l};
                        sort(ij, ij + 2);
                        sort(ik, ik + 2);
                        sort(il, il + 2);
                        saturate_z.push_back(&(zzpm[get_hash({ij[1], ij[0], k, l})]));
                        saturate_z.push_back(&(zzpm[get_hash({ij[1], ij[0], l, k})]));
                        saturate_z.push_back(&(zzpm[get_hash({ik[1], ik[0], j, l})]));
                        saturate_z.push_back(&(zzpm[get_hash({ik[1], ik[0], l, j})]));
                        saturate_z.push_back(&(zzpm[get_hash({il[1], il[0], k, j})]));
                        saturate_z.push_back(&(zzpm[get_hash({il[1], il[0], j, k})]));
                    }
                }

            }
        }

    }

    for (auto m : saturate_z) {
        (*m) = complex0;
    }

}

void Density_matrix::store_br(double **br_zz, double **br_pm, double **br_mp, const Crystal &crystal) {
	/* Temporarily, using br_zz as 2 spin
	 * br_pm as three spin
	 * br_mp as four spin
	 */
    for (int i = 0; i < number_of_spins; i++) {
        for (auto j : crystal.nearest_neighbours[i]) {
            // two spin terms
            //br_zz[j][i] += abs(zz[get_hash({j, i})]);
            br_zz[j][i] += abs(pm[get_hash({i, j})]);
            //br_mp[j][i] += abs(pm[get_hash({j, i})]);

            // three spin
            for (int k = 0; k < crystal.number_of_nuclei; k++) {
                int ijk[3] = {i, j, k};
                sort(ijk, ijk + 3);
                int hash;
                //hash = get_hash({ijk[2], ijk[1], ijk[0]});
                //if (zzz.find(hash) != zzz.end()) br_zz[j][i] += abs(zzz[hash]);
                hash = get_hash({k, i, j});
                if (zpm.find(hash) != zpm.end()) br_pm[j][i] += abs(zpm[hash]);
                //hash = get_hash({k, j, i});
                //if (zpm.find(hash) != zpm.end()) br_mp[j][i] += abs(zpm[hash]);

                // four spin
                for (int l = 0; l < crystal.number_of_nuclei; l++) {
                    // zzpm
                    // ppmm

                    // zzpm
                    //hash = get_hash({j, i, k, l});
                    //if (zzpm.find(hash) != zzpm.end()) br_zz[j][i] += abs(zzpm[hash]);
                    //hash = get_hash({j, i, l, k});
                    //if (zzpm.find(hash) != zzpm.end()) br_zz[j][i] += abs(zzpm[hash]);

                    //int kl[2] = {k, l};
                    //sort(kl, kl + 2);
                    //hash = get_hash({kl[1], kl[0], i, j});
                    //if (zzpm.find(hash) != zzpm.end()) br_pm[j][i] += abs(zzpm[hash]);
                    //hash = get_hash({kl[1], kl[0], j, i});
                    //if (zzpm.find(hash) != zzpm.end()) br_mp[j][i] += abs(zzpm[hash]);

                    // ppmm
                    // we only want ZQ. So ppmm[i, .., .., j], ppmm[j, .., .., i]

                    int ik[2] = {i, k};
                    sort(ik, ik + 2);
                    int il[2] = {i, l};
                    sort(il, il + 2);
                    int jk[2] = {j, k};
                    sort(jk, jk + 2);
                    int jl[2] = {j, l};
                    sort(jl, jl + 2);

                    // I+, J-
                    hash = get_hash({ik[1], ik[0], jl[1], jl[0]});
                    if (ppmm.find(hash) != ppmm.end()) br_mp[j][i] += abs(ppmm[hash]);
                    hash = get_hash({il[1], il[0], jk[1], jk[0]});
                    if (ppmm.find(hash) != ppmm.end()) br_mp[j][i] += abs(ppmm[hash]);

                    // I-, J+
                    //hash = get_hash({jk[1], jk[0], il[1], il[0]});
                    //if (ppmm.find(hash) != ppmm.end()) br_mp[j][i] += abs(ppmm[hash]);
                    //hash = get_hash({jl[1], jl[0], ik[1], ik[0]});
                    //if (ppmm.find(hash) != ppmm.end()) br_mp[j][i] += abs(ppmm[hash]);
                }
            }
        }
    }
}
