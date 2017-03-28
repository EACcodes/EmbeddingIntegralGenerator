/*
 * @file TigerStructs.h
 * @brief Sets up tiger namespace for generalization of basis set and integral data
 *
 * Author: Francis Ricci
 * Date created: 1/26/15
 *
 */


#ifndef TIGER_STRUCTS
#define TIGER_STRUCTS

#include <armadillo>

namespace tiger
{
enum basis_type {full, minimal, aux};
enum return_val {SUCCESS = 0, INPUT_FAILURE = -1, MO_PARSE_FAILURE = -2, BASIS_FAILURE = -3,
	ORBITAL_FAILURE = -4, CHOLESKY_FAILURE = -5, DF_FAILURE = -6, TIGER_FAILURE = -7,
    EMB_INT_FAILURE = -8};

typedef struct {
	double x;
	double y;
	double z;
} coords_t;

typedef struct {
	coords_t coords;
	int Z;
} atom_t;

typedef struct {
    size_t l;
    size_t m;
    size_t n;
} ang_mom_comp_t;

typedef struct {
    coords_t center;
    size_t ind_center;
    size_t ang_mom;
    size_t nbf;
    bool use_sph;
    std::vector<double> coeff;
    std::vector<double> coeff_norm;
    std::vector<double> exp;
    size_t first_ind;
    std::vector<double> relnorms;
    std::vector<ang_mom_comp_t> ang_mom_comp;
} shell_t;

typedef struct {
    size_t ang_mom;
    size_t n_prim;
    size_t n_contr;
    arma::mat contr_mat;
    std::vector<double> exp;
    size_t first_ind;
} molcas_shell_t;


}
#endif
