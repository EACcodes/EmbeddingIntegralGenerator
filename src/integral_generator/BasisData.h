/*
 * @file BasisData.h
 *
 *  Created on: Jan 15, 2015
 *      Author: Francis Ricci
 */

#ifndef TIGER_CPP_BASIS_H_
#define TIGER_CPP_BASIS_H_

#include "ConvertKeys.h"
#include "TigerStructs.h"

#include <armadillo>
#include <vector>
#include <memory>
#include <list>

/**
 * @class BasisData
 * @brief An abstract class, presenting an interface for basis set and integral data
 */

class BasisData {
public:
	//copy constructor
	virtual BasisData* duplicate() = 0;
	virtual BasisData* decontract() = 0;  // construct decontracted basis set

	//basis data
	virtual std::vector<tiger::coords_t> get_coordinates() = 0;
	virtual std::vector<tiger::atom_t> get_atoms() = 0;
	virtual std::vector<int> get_nFuncInShell() = 0;
	virtual int get_nBas() = 0;
	virtual int get_nShell() = 0;
	virtual int get_Ncart() = 0;       //
	virtual int get_nAtoms() = 0;
	virtual std::vector<int> get_basis_map() = 0; //maps basis functions to atom number
	virtual arma::vec eval_func(double x, double y, double z) = 0; // evaluates the basis functions at point (x,y,z)
	virtual arma::vec eval_func(size_t shell_ind, double x, double y, double z) = 0; // evaluates the functions of shell shell_ind at point (x,y,z)
	virtual void print_basis() = 0; // print basis functions

	// Shell information
	virtual void print_shell(size_t ind) = 0; // print shell with index idn
	virtual std::vector<tiger::shell_t> get_shells() = 0; // Get Shells (coefficients, exponents, coords, etc.)
	virtual arma::mat get_trans_mat() = 0; // Get transformation matrix between contracted and uncontracted basis functions
	virtual bool is_spherical(size_t ind) = 0; // Are spherical harmonics in use?
	virtual arma::mat get_trans_cart_sph() = 0; // Get transformation matrix from cartesian to spherical
	virtual arma::mat get_trans_cart_sph(size_t ind) = 0; // Get transformation matrix from cartesian to spherical for shell index ind
	virtual size_t get_shell_center_ind(size_t ind) = 0;

	virtual ~BasisData();

	ConvertKeys *_globals;

protected:
	tiger::basis_type _type;
};

/**
 * @class EmbeddingBasis
 * @brief Uses own implementation of the BasisData interface for basis set and integral data
 */

class EmbeddingBasis : public BasisData{
public:
	EmbeddingBasis(ConvertKeys& global_vars);
	~EmbeddingBasis();

	//copy constructor (for thread-safe integrals)
	BasisData* duplicate();
	BasisData* decontract();
	EmbeddingBasis(EmbeddingBasis& bas,bool decontracted);

	//basis data
	std::vector<tiger::coords_t> get_coordinates();
	std::vector<tiger::atom_t> get_atoms();
	std::vector<int> get_nFuncInShell();
	int get_nBas();
	int get_nShell();
	int get_Ncart();
	int get_Nsph();
	int get_nAtoms();
	std::vector<int> get_basis_map();
	arma::vec eval_func(double x, double y, double z);
	arma::vec eval_func(size_t shell_ind, double x, double y, double z);
	void print_basis();

	// Shell information
	void print_shell(size_t ind);
	bool is_spherical(size_t ind);
	arma::mat get_trans_cart_sph(size_t ind);
	arma::mat get_trans_cart_sph();
	std::vector<tiger::shell_t> get_shells();
	arma::mat get_trans_mat();
	size_t get_shell_center_ind(size_t ind);

protected:
	std::vector<tiger::shell_t> _shells;
	std::vector<tiger::atom_t> _atoms;

	arma::mat _trans; // transformation matrix between decontracted and original basis (nr primitives x nr basis functions)
	const double _ang2bohr = 1./0.529177249;


	std::vector<tiger::atom_t> read_structure(string file_name); // read geometric data from xyz file
	void create_basis();
	std::list<string> read_basis_file(string file_name);  // read basis set file into  list of lines
														  // delete comments, empty lines and convert to lower case
	std::vector<tiger::molcas_shell_t> read_shells(double nuc_charge, std::list<string> &lines);
	int dfact(int number);
	int fact(int number);
	int bin_coeff(int n, int k);
	arma::mat trans_cart_sph(int l);
	arma::mat trans_cart_sph_table(int l);
	arma::mat trans_cart_sph_old(int l);
	int nuc_chg(std::string token);
};

#endif /* TIGER_CPP_BASIS_H_ */
