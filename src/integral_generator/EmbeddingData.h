/*
 * @file EmbeddingData.h
 *
 *      Author: Caroline M. Krauter
 */

#ifndef TIGER_CPP_EMBEDDING_H_
#define TIGER_CPP_EMBEDDING_H_

#include "integral_generator/ConvertKeys.h"
#include "integral_generator/BasisData.h"

#include <armadillo>
#include <vector>

/**
 * @class EmbeddingData
 * @brief Contains the implementation of embedding data 
 */

typedef struct {
    double x;
    double y;
    double z;
    double factor;
} kpoint_t;

class EmbeddingData {
public:
    // constructor
    EmbeddingData(ConvertKeys& global_vars, BasisData *bas);
    ~EmbeddingData();

    // Compute the one-electron integrals over the embedding potential
    void compute_integrals();

private:
    //generate embedding data
    arma::mat compute_emb_one_ints_real();
    double compute_cell_volume(std::vector<std::vector<double>>& lattice_vecs);
    int compute_nr_grid_pts(std::vector<int>& grid_dimension);
    std::vector<double> read_embedding_potential(int nr_grid_points);
    std::vector<std::vector<double>> convert_to_bohr(std::vector<std::vector<double>>& vecs);
    std::vector<double> cross_product(std::vector<double> v1, std::vector<double> v2, double factor);
    std::vector<complex<double>> int_gauss_pw(double x, double l, double alpha);
    std::vector<kpoint_t> construct_kpoints(std::vector<std::vector<double>>& reciprocal_vecs_incr);
    std::vector<std::vector<double>> construct_grid_real(std::vector<std::vector<double>>& incr_lattice_vecs);
    arma::mat compute_real_space_part();
    arma::mat compute_reciprocal_space_part();

    //data
    BasisData *_basis;
    ConvertKeys *_globals;
    std::string _method;

    // grid information
    std::vector<int> _grid_dimension;
    std::vector<std::vector<double>> _shift_vec;
    std::vector<std::vector<double>> _lattice_vecs;
    int _nr_grid_points;
    double _cell_volume;

    // embedding potential
    std::vector<double> _embedding_potential;

    // Use this conversion factor for the moment for better comparability with Molcas
    // TODO: Eventually change this to the one used in the rest of the code
    //
    // from embedding implementation in Molcas. Originally from (confirmed):
    // *     Conversion factor angstrom to bohr from the IUPAC
    // *     publication
    // *     .529177249(24) angstrom / bohr
    // *     "Quantities, Units and Symbols in Physical Chemistry"
    // *     I. Mills, T. Cvitas, K. Homann, N. Kallay and
    // *     K. Kuchitsu, Blackwell Scientific Publications,
    // *     Oxford, 1988.
    const double _ang2bohr = 1./0.529177249;

    // Cut-off value for sum of exponents; if larger do integral on reciprocal grid (as in MOLCAS version by Kuang)
    const double _alpha_cut = 20.;

    // Cut-off value (as in MOLCAS version by Kuang)
    const double _exp_thresh = 20.;

    // print embedding potentials
    void print_tigerci(arma::mat ints, std::string filename);
    void print_molpro(arma::mat ints, std::string filename);
    void print_gamess(arma::mat ints, std::string filename);
    void print_molcas(arma::mat ints, std::string filename);
    void print_pyscf(arma::mat ints, std::string filename);

    // sort to Molpro order
    arma::mat sort_molpro(arma::mat ints);
    // maps for resorting functions
    // TigerCI: AOs are sorted with respect to -l,-l+1,...,0,1,...,l-1,l for each l>1, 
    //           for l<=1: s,px,py,pz
    // Molpro:
    //  '1s',
    //  '2px','2py','2pz'
    //  '3d0','3d2-','3d1+','3d2+','3d1-'
    //  '4f1+','4f1-','4f0','4f3+','4f2-'
    //  '4f3-','4f2+'
    //  '5g0','5g2-','5g1+','5g4+','5g1-','5g2+'
    //  '5g4-','5g3+','5g3-'
    //  '6h1+','6h1-','6h2+','6h3+','6h4-','6h3-','6h4+','6h5-','6h0','6h5+','6h2-'
    //  '7i6+','7i2-','7i5+','7i4+','7i5-','7i2+','7i6-','7i3+','7i4-','7i0','7i3-','7i1-','7i1+'
    // Maps for different angular momenta; numbers give Molpro-index of this TigerCI-component
    int _map_s[1] = {0};
    int _map_p[3] = {0,1,2};
    int _map_d[5] = {1,4,0,2,3};
    int _map_f[7] = {5,4,1,2,0,6,3};
    int _map_g[9] = {6,8,1,4,0,2,5,7,3};
    int _map_h[11] = {7,4,5,10,1,8,0,2,3,6,9};
    int _map_i[13] = {6,4,8,10,1,11,9,12,5,7,3,2,0};

    arma::mat sort_gamess(arma::mat ints);
    // maps for resorting functions
    // TigerCI: Use Cartesian functions
    //           s: s
    //           p: x, y, z
    //           d: xx, xy, xz, yy, yz, zz
    //           f: xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz
    //
    // GAMESS (l_max=3!!!):
    // TDOD: Throw error before calculation if GAMESS requested together with l>3
    //   S
    //   X, Y, Z
    //   XX , YY ,  ZZ ,  XY ,  XZ ,  YZ
    //   XXX , YYY , ZZZ , XXY , XXZ , YYX , YYZ , ZZX , ZZY , XYZ
    // Maps for different angular momenta; numbers give GAMESS-index of this TigerCI-component
    int _map_gamess_s[1] = {0};
    int _map_gamess_p[3] = {0,1,2};
    int _map_gamess_d[6] = {0,3,4,1,5,2};
    int _map_gamess_f[10] = {0,3,4,5,9,7,1,6,8,2};


};

#endif /* TIGER_CPP_EMBEDDING_H_ */
