/*
 * @file EmbeddingData.cpp
 *
 *      Author: Caroline M. Krauter
 */

#include "integral_generator/EmbeddingData.h"
#include "integral_generator/Timer.h"
#include "integral_generator/BasisData.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <math.h>
#include <complex>

#include <fftw3.h>

EmbeddingData::~EmbeddingData(){}

EmbeddingData::EmbeddingData(ConvertKeys& global_vars, BasisData *bas) {
    _globals = &global_vars;
    _basis = bas;

    Timer embedding_info;

    cout << "***********************" << endl;
    cout << "*                      " << endl;
    cout << "* Embedding Integrals  " << endl;
    cout << "*                      " << endl;
    cout << "***********************" << endl;
    cout << endl;

    // get information about cell from input file
    std::vector<std::vector<int>> vecs = _globals->get_global_int_vecs("grid dimensions");
    _grid_dimension = vecs[0];
    _shift_vec = _globals->get_global_double_vecs("shift vector");
    _lattice_vecs = _globals->get_global_double_vecs("lattice vectors");

    // convert to bohr
    _shift_vec = convert_to_bohr(_shift_vec);
    _lattice_vecs = convert_to_bohr(_lattice_vecs);

    // calculate number of grid points and cell volume
    _nr_grid_points = compute_nr_grid_pts(_grid_dimension);
    _cell_volume = compute_cell_volume(_lattice_vecs);

    _embedding_potential = read_embedding_potential(_nr_grid_points);
    _method = _globals->get_global_string("output format");

    embedding_info.print_elapsed_time("  Read grid data and embedding potential");
}

void EmbeddingData::compute_integrals() {

    Timer embedding_integrals;

    // Initialize matrix for one electron integrals
    arma::mat ints;

    // Do integration either in real space or ...
    if (_globals->get_global_bool("embedding ints real")) {
        cout << endl << "* Real space integration scheme  " << endl;
        ints = compute_emb_one_ints_real();

        embedding_integrals.print_elapsed_time("  Calculate embedding integrals");
    }
    // ... with the reciprocal space scheme
    else {
    	cout << endl << "* Reciprocal space integration scheme" << endl;
        // real-space part
        ints = compute_real_space_part();
        embedding_integrals.print_elapsed_time("  Calculate embedding integrals, real space part");

        // reciprocal space part
        ints += compute_reciprocal_space_part();

        embedding_integrals.print_elapsed_time("  Calculate embedding integrals, reciprocal space part");
    }

    // Do contractions and transformation to spherical basis functions as needed and print results
    if (_method == "TIGERCI") {
    	// contract
    	arma::mat trans_contr = _basis->get_trans_mat();
    	ints = trans_contr.t() * ints * trans_contr ;
    	// transform to spherical
    	arma::mat trans_sph = _basis->get_trans_cart_sph();
    	ints = trans_sph * ints * trans_sph.t();

    	print_tigerci(ints,"emb_ints.tigerci");
    }
    else if (_method == "MOLPRO") {
    	// contract
    	arma::mat trans_contr = _basis->get_trans_mat();
    	ints = trans_contr.t() * ints * trans_contr ;
    	// transform to spherical
    	arma::mat trans_sph = _basis->get_trans_cart_sph();
    	ints = trans_sph * ints * trans_sph.t();

        print_molpro(sort_molpro(ints),"emb_ints.molpro");
    }
    else if (_method == "PYSCF") {
        // contract
        arma::mat trans_contr = _basis->get_trans_mat();
        ints = trans_contr.t() * ints * trans_contr ;
        // transform to spherical
        arma::mat trans_sph = _basis->get_trans_cart_sph();
        ints = trans_sph * ints * trans_sph.t();

        print_pyscf(ints,"emb_ints.pyscf");
    }
    else if (_method == "GAMESS") {
    	// contract
    	arma::mat trans_contr = _basis->get_trans_mat();
    	ints = trans_contr.t() * ints * trans_contr ;

    	print_gamess(sort_gamess(ints),"emb_ints.gamess");
    }
    else if (_method == "NWCHEM") {
    	// contract
    	arma::mat trans_contr = _basis->get_trans_mat();
    	ints = trans_contr.t() * ints * trans_contr ;

    	cout << "NWCHEM does not require resorting of the integrals, but a specific output format has not yet been implemented. Printing TigerCI to save the result" << endl;
        print_tigerci(ints,"emb_ints.backup.tigerci");
    }
    else if (_method == "MOLCAS") {
    	print_molcas(ints,"emb_ints.molcas");
    }
    else if (_method == "QCHEM") {
	//same as tigerci
        // contract
        arma::mat trans_contr = _basis->get_trans_mat();
        ints = trans_contr.t() * ints * trans_contr ;
        // transform to spherical
        arma::mat trans_sph = _basis->get_trans_cart_sph();
        ints = trans_sph * ints * trans_sph.t();

        print_tigerci(ints,"emb_ints.qchem");
	_basis->print_basis();
	_basis->print_basis_gaussian94_format();

    }
    else {
        cout << "This output format is not recognized and no output will be generated. ";
        cout << "This should have been detected when reading the input. Please let a developer know.";
    }
    cout << endl;
    embedding_integrals.print_elapsed_time("  Calculate embedding integrals, sorting and printing");
}

arma::mat EmbeddingData::compute_emb_one_ints_real() {

    Timer embedding_integrals;

    // determine incremental lattice vectors in bohr
    std::vector<std::vector<double>> incr_lattice_vecs (3, std::vector<double>(3));
    for (size_t i=0; i<3; i++){
        for (size_t j =0; j<3; j++){
            incr_lattice_vecs[i][j] = _lattice_vecs[i][j]/_grid_dimension[i];
        }
    }

    // determine volume element for integration
    double volume_element = _cell_volume / _nr_grid_points;

    // Debug output
    cout << "nr_grid_points " << _nr_grid_points << endl;
    cout << std::setprecision (15) << "cell_volume " << _cell_volume << endl;
    cout << std::setprecision (15) << "volume_elemet " << volume_element << endl;
    cout << "incremental lattice vectors in bohr: " << endl;
    for (size_t i=0; i<incr_lattice_vecs.size();i++){
        for (size_t j =0; j<incr_lattice_vecs[i].size();j++){
            cout << incr_lattice_vecs[i][j] << " ";
        } 
        cout << endl;
    }

    BasisData *decontr_basis = _basis->decontract();
    std::vector<tiger::shell_t> decontr_shells  = decontr_basis->get_shells();

    // initiate matrix for integrals over primitive basis functions
    int nbas_prim = decontr_basis->get_nBas();
    int nshell_prim = decontr_basis->get_nShell();
    arma::mat ints_prim;
    ints_prim.zeros(nbas_prim,nbas_prim);

    // loop over grid points
    size_t i = 0; // index to determine the index of a grid point in V_emb
    for(size_t i_z=0; i_z<_grid_dimension[2]; i_z++) {
        for (size_t i_y=0; i_y<_grid_dimension[1]; i_y++) {
            for (size_t i_x=0; i_x<_grid_dimension[0]; i_x++, i++) {

            	// determine coordinates of grid point
                double x = (_shift_vec[0][0] + i_x * incr_lattice_vecs[0][0] + i_y * incr_lattice_vecs[1][0] + i_z * incr_lattice_vecs[2][0]);
                double y = (_shift_vec[0][1] + i_x * incr_lattice_vecs[0][1] + i_y * incr_lattice_vecs[1][1] + i_z * incr_lattice_vecs[2][1]);
                double z = (_shift_vec[0][2] + i_x * incr_lattice_vecs[0][2] + i_y * incr_lattice_vecs[1][2] + i_z * incr_lattice_vecs[2][2]);

            	// evaluate basis functions at this point
                arma::vec funcs = decontr_basis->eval_func(x, y, z);

                // add products (\phi_a * V_emb * \phi_b * volume_element) to the matrix of integrals
                // use symmetry of integrals: int(a,b)=int(b,a)
                for (size_t a=0; a<nbas_prim; a++) {
                    ints_prim(a,a) += funcs(a) * funcs(a) * volume_element * _embedding_potential[i];
                    for (size_t b=0; b<a; b++) {
                        ints_prim(a,b) += funcs(a) * funcs(b) * volume_element * _embedding_potential[i];
                    }
                }
            }
        }
    }

    // add remaining elements based on symmetry: int(b,a)=int(a,b)
    for (size_t a=0; a < nbas_prim; a++) {
        for (size_t b =0; b < a; b++) {
            ints_prim(b,a)=ints_prim(a,b);
        }
    }

    return ints_prim;
}

int EmbeddingData::compute_nr_grid_pts(std::vector<int>& grid_dimension) {
    return grid_dimension[0]*grid_dimension[1]*grid_dimension[2];
}

double EmbeddingData::compute_cell_volume(std::vector<std::vector<double>>& lattice_vecs) {

    // volume of cell spanned by lattice vectors a,b,c is:
    // |a * (b x c)|

    // calculate cross product beween b and c
    std::vector<double> cross_prod = cross_product(lattice_vecs[1],lattice_vecs[2],1.0);

    // Calculate volume 
    return std::abs(cross_prod[0] * lattice_vecs[0][0] + cross_prod[1] * lattice_vecs[0][1] + cross_prod[2] * lattice_vecs[0][2]);
}

std::vector<double> EmbeddingData::read_embedding_potential(int nr_grid_points) {

    string file_name = _globals->get_global_string("embedding potential file");
    ifstream inputf(file_name);
    if (!inputf){
    	throw runtime_error("Problem opening embedding potential file.");
    }

    std::vector<double> embedding_potential;
    embedding_potential.reserve(nr_grid_points);

    double number;
    while (inputf >> number) {
        embedding_potential.push_back(number);
    }
    // check whether the correct number of values was read
    if (embedding_potential.size() != nr_grid_points) {
        std::ostringstream oss;
        oss << "Problem reading embedding potential: Number of values read (= "<< embedding_potential.size() <<
         		") is not the same as number of grid points (=" << nr_grid_points << ")!\n";
        throw std::runtime_error(oss.str());
    }
    return embedding_potential;
}

std::vector<std::vector<double>> EmbeddingData::convert_to_bohr(std::vector<std::vector<double>>& vecs) {
    std::vector<std::vector<double>> vecs_conv = vecs;
    for (size_t i=0; i<vecs.size(); i++){
        for (size_t j =0; j<vecs[i].size(); j++){
            vecs_conv[i][j] = vecs[i][j] * _ang2bohr;
        }
    }
    return vecs_conv;
}   

std::vector<double> EmbeddingData::cross_product(std::vector<double> v1, std::vector<double> v2, double factor) {
    std::vector<double> cross_prod;
    cross_prod.push_back((v1[1] * v2[2] - v1[2] * v2[1])*factor);
    cross_prod.push_back((v1[2] * v2[0] - v1[0] * v2[2])*factor);
    cross_prod.push_back((v1[0] * v2[1] - v1[1] * v2[0])*factor);
    return cross_prod;
}

std::vector<complex<double>> EmbeddingData::int_gauss_pw(double x, double l, double alpha) {
    // calculate  \int e^ikx * x^m * exp^(-alpha*x^2) for m=l,l-1,...,0
	// Yu, K.; Libisch, F. & Carter, E. A. The Journal of Chemical Physics, 2015, 143
	//
	// f(x,0,a) = sqrt(pi/a) * exp(-x**2/4a) #remark: Here, this is treated as a prefactor and calculated elsewhere
	// f(x,1,a) = ik/2a * f(x,0,a)
	// f(x,l,a) = (l-1)/2a * f(x,l-2,a) + ik/2a * f(x,l-1,a)

    std::vector<complex<double>> result;
    
    for (size_t m = 0 ; m <= l; m++) {
        if (m==0) result.push_back(complex<double>(1.0,0.0));
        else if (m==1) result.push_back(complex<double>(0,1) * x * 0.5 / alpha * result[0]);
        else result.push_back( ( (m-1) * 0.5 * result[m-2] +  complex<double>(0,1) * x * 0.5 * result[m-1] ) / alpha );
    }

    return result;

}

std::vector<kpoint_t> EmbeddingData::construct_kpoints(std::vector<std::vector<double>>& reciprocal_vecs_incr) {

    std::vector<kpoint_t> kpoints;
    for(size_t i_z=0; i_z<_grid_dimension[2]; i_z++) {
        int kz = i_z;
        if (kz > _grid_dimension[2]/2) {kz = kz - _grid_dimension[2];}           

        for (size_t i_y=0; i_y<_grid_dimension[1]; i_y++) {
            int ky = i_y;
            if (ky > _grid_dimension[1]/2) {ky = ky - _grid_dimension[1];}                                                                                     
            for (size_t i_x=0; i_x<(_grid_dimension[0]/2+1); i_x++) {
                int kx = i_x;
                kpoint_t kpoint;

                // determine prefactor
                if (i_x == 0 || _grid_dimension[0] - i_x == i_x) {
                    kpoint.factor = 1.0;
                }
                else {                                                                                                                                                                                                          
                    kpoint.factor = 2.0;
                }

                // determine coordinates of reciprocal grid point                                                                                           
                kpoint.x = (kx * reciprocal_vecs_incr[0][0] + ky * reciprocal_vecs_incr[1][0] + kz * reciprocal_vecs_incr[2][0]);   
                kpoint.y = (kx * reciprocal_vecs_incr[0][1] + ky * reciprocal_vecs_incr[1][1] + kz * reciprocal_vecs_incr[2][1]);   
                kpoint.z = (kx * reciprocal_vecs_incr[0][2] + ky * reciprocal_vecs_incr[1][2] + kz * reciprocal_vecs_incr[2][2]);   

                kpoints.push_back(kpoint);
            }
        }
    }
    return kpoints;
}

std::vector<std::vector<double>> EmbeddingData::construct_grid_real(std::vector<std::vector<double>>& incr_lattice_vecs) {
    std::vector<std::vector<double>> points;
    for(size_t i_z=0; i_z<_grid_dimension[2]; i_z++) {
        for (size_t i_y=0; i_y<_grid_dimension[1]; i_y++) {
            for (size_t i_x=0; i_x<_grid_dimension[0]; i_x++) {
                std::vector<double> x(3);
                for (size_t k=0; k<3; k++) {
                    x[k] = _shift_vec[0][k] + i_x * incr_lattice_vecs[0][k] + i_y * incr_lattice_vecs[1][k] + i_z * incr_lattice_vecs[2][k];
                }
                points.push_back(x);
            }
        }
    }
    return points;
}

arma::mat EmbeddingData::compute_real_space_part() {

    // determine volume element for integration
    double volume_element = _cell_volume / _nr_grid_points;

    // decontract basis and get shells 
    BasisData *decontr_basis = _basis->decontract();
    std::vector<tiger::shell_t> decontr_shells  = decontr_basis->get_shells();
    // initiate matrix for integrals over primitive basis functions
    int nbas_prim = decontr_basis->get_nBas();
    int nshell_prim = decontr_basis->get_nShell();
    arma::mat ints_prim;
    ints_prim.zeros(nbas_prim,nbas_prim);

    // determine incremental lattice vectors in bohr
    std::vector<std::vector<double>> incr_lattice_vecs (3, std::vector<double>(3));
    for (size_t i=0; i<3; i++){
        for (size_t j =0; j<3; j++){
            incr_lattice_vecs[i][j] = _lattice_vecs[i][j]/_grid_dimension[i];
        }
    }

    //Construct incremental reciprocal vectors
    std::vector<std::vector<double>> reciprocal_vecs_incr;
    reciprocal_vecs_incr.push_back(cross_product(_lattice_vecs[1],_lattice_vecs[2],2.0*M_PI/_cell_volume));
    reciprocal_vecs_incr.push_back(cross_product(_lattice_vecs[2],_lattice_vecs[0],2.0*M_PI/_cell_volume));
    reciprocal_vecs_incr.push_back(cross_product(_lattice_vecs[0],_lattice_vecs[1],2.0*M_PI/_cell_volume));

    // Construct real space grid
    std::vector<std::vector<double>> points = construct_grid_real(incr_lattice_vecs);

    // Calculate bounds used to determine whether box large enough for integrals over shell pairs
    std::vector<double> bounds;
    for (size_t i=0; i<3; i++){
        bounds.push_back( M_PI * M_PI / (
            reciprocal_vecs_incr[i][0] * reciprocal_vecs_incr[i][0] +
            reciprocal_vecs_incr[i][1] * reciprocal_vecs_incr[i][1] +
            reciprocal_vecs_incr[i][2] * reciprocal_vecs_incr[i][2] ));
    }

    /////////////////////////////
    //Do real space integration//
    /////////////////////////////

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nshell_prim,bounds,decontr_basis,points,decontr_shells,reciprocal_vecs_incr,ints_prim,volume_element,cout)
#endif
    // loop over shells
    for (size_t a=0; a<nshell_prim; a++) {
    	double alpha_a = decontr_shells[a].exp[0];
        for (size_t b=0; b<=a; b++) {
            double alpha_b = decontr_shells[b].exp[0];

            double exp1 = (pow((decontr_shells[a].center.x-decontr_shells[b].center.x),2) +
               pow((decontr_shells[a].center.y-decontr_shells[b].center.y),2) +
               pow((decontr_shells[a].center.z-decontr_shells[b].center.z),2)) *
               alpha_a * alpha_b / (alpha_a + alpha_b);

            if ((decontr_shells[a].ind_center == decontr_shells[b].ind_center) && (alpha_a + alpha_b > _alpha_cut)) {
                continue; // will be done in reciprocal space
            }
            else if ( exp1 > _exp_thresh ) {
                continue; // skip this shell pair if (r_a-r_b)**2 * alpha_a * alpha_b / (alpha_a + alpha_b) > _exp_thresh
            }
            else {

                // Center of combined Gaussian
                std::vector<double> RI;
                RI.push_back((decontr_shells[a].center.x * alpha_a + decontr_shells[b].center.x * alpha_b) / (alpha_a + alpha_b));
                RI.push_back((decontr_shells[a].center.y * alpha_a + decontr_shells[b].center.y * alpha_b) / (alpha_a + alpha_b));
                RI.push_back((decontr_shells[a].center.z * alpha_a + decontr_shells[b].center.z * alpha_b) / (alpha_a + alpha_b));

                // Padding necessary in any direction? -> set ranges to 1 for that direction
                std::vector<int> ranges(3,0);
                for (size_t k = 0 ; k<3; k++) {
                    if (( alpha_a + alpha_b ) * bounds[k] < _exp_thresh) {
                        ranges[k]=1;
                    }
                }

                arma::mat shell_ints(decontr_shells[a].ang_mom_comp.size(),decontr_shells[b].ang_mom_comp.size());
                shell_ints.zeros();

                // loop over grid points
                for(size_t i=0; i<_nr_grid_points; i++) {

                    // Periodic Boundary condition shifts (shift so that center of combined Gaussian in the center)
                    std::vector<double> dr(3);
                    for (size_t k=0; k<3; k++) {
                        dr[k] = points[i][k] - RI[k];
                    }
                    std::vector<double> ds(3);
                    for (size_t k=0; k<3; k++) {
                        ds[k] = (reciprocal_vecs_incr[k][0] * dr[0] + reciprocal_vecs_incr[k][1] * dr[1] + reciprocal_vecs_incr[k][2] * dr[2]) / 2.0 / M_PI ;
                        ds[k]-= floor(ds[k]+0.5);
                    }


                    // if padded, loop over adjacent cells
                    std::vector<double> dsp(3);
                    for (int ii_z = -ranges[2]; ii_z < ranges[2]+1; ii_z++) {
                        dsp[2] = ds[2] + ii_z;
                        for (int ii_y = -ranges[1]; ii_y < ranges[1]+1; ii_y++) {
                            dsp[1] = ds[1] + ii_y;
                            for (int ii_x = -ranges[0]; ii_x < ranges[0]+1; ii_x++) {
                                dsp[0] = ds[0] + ii_x;

                                // position of point wrt center of combined Gaussian
                                for (size_t k=0; k<3; k++) {
                                    dr[k] = dsp[0] * _lattice_vecs[0][k] + dsp[1] * _lattice_vecs[1][k] + dsp[2] * _lattice_vecs[2][k] ;

                                }

                                double exp2 = (dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2])*(alpha_a + alpha_b);
                                if (exp2 > _exp_thresh) {
                                    continue;
                                } // skip if exponent of combined Gaussian at this point very large -> value approx. 0

                                // points where functions need to be evaluated
                                std::vector<double> xA(3);
                                std::vector<double> xB(3);
                                for (size_t k=0; k<3; k++) {
                                    xA[k] = dr[k] + RI[k];
                                    xB[k] = dr[k] + RI[k];
                                }

                                // Coordinates relative to center of Gaussians
                                xA[0] = xA[0] - decontr_shells[a].center.x;
                                xA[1] = xA[1] - decontr_shells[a].center.y;
                                xA[2] = xA[2] - decontr_shells[a].center.z;

                                xB[0] = xB[0] - decontr_shells[b].center.x;
                                xB[1] = xB[1] - decontr_shells[b].center.y;
                                xB[2] = xB[2] - decontr_shells[b].center.z;

                                // angular independent part of integral
                                double radial = exp(-(exp1+exp2)) *  _embedding_potential[i];

                                // multiply with angular parts
                                for (size_t ai = 0; ai < decontr_shells[a].ang_mom_comp.size() ; ai++) {
                                    double ang_a =
                                    		pow(xA[0],decontr_shells[a].ang_mom_comp[ai].l) *
                                    		pow(xA[1],decontr_shells[a].ang_mom_comp[ai].m) *
                                    		pow(xA[2],decontr_shells[a].ang_mom_comp[ai].n);
                                    for (size_t bi = 0; bi < decontr_shells[b].ang_mom_comp.size() ; bi++) {
                                    	double ang_b =
                                        		pow(xB[0],decontr_shells[b].ang_mom_comp[bi].l) *
                                        		pow(xB[1],decontr_shells[b].ang_mom_comp[bi].m) *
                                        		pow(xB[2],decontr_shells[b].ang_mom_comp[bi].n);
                                        shell_ints(ai,bi) += radial * ang_a * ang_b;
                                    }
                                }

                            }
                        }
                    }
                }
                for (size_t ai = 0; ai < decontr_shells[a].nbf ; ai++) {
                	for (size_t bi = 0; bi < decontr_shells[b].nbf ; bi++) {
                		ints_prim(decontr_shells[a].first_ind+ai , decontr_shells[b].first_ind+bi)= shell_ints(ai,bi)
                				*decontr_shells[a].coeff_norm[0] * decontr_shells[b].coeff_norm[0]
                				*decontr_shells[a].relnorms[ai] *decontr_shells[b].relnorms[bi] * volume_element;
                	}
                }
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nbas_prim,ints_prim)
#endif

    // add remaining elements based on symmetry: int(b,a)=int(a,b)
    for (size_t a=0; a < nbas_prim; a++) {
        for (size_t b =0; b < a; b++) {
            ints_prim(b,a)=ints_prim(a,b);
        }
    }

    return ints_prim;

}


arma::mat EmbeddingData::compute_reciprocal_space_part() {

    // decontract basis and get shells 
    BasisData *decontr_basis = _basis->decontract();
    std::vector<tiger::shell_t> decontr_shells  = decontr_basis->get_shells();
    // initiate matrix for integrals over primitive basis functions
    int nbas_prim = decontr_basis->get_nBas();
    int nshell_prim = decontr_basis->get_nShell();
    arma::mat ints_prim;
    ints_prim.zeros(nbas_prim,nbas_prim);

    // Prepare FFT 
    double* pot = &_embedding_potential[0];
    fftw_complex *pot_fft;
    // nr grid points of transformed grid:
    size_t nr_grid_points_trans = ((_grid_dimension[0]/2 +1) * _grid_dimension[1] * _grid_dimension[2]);

    //Do FFT
    pot_fft = (fftw_complex*) fftw_malloc(nr_grid_points_trans * sizeof(fftw_complex));
    fftw_plan plan;
    plan = fftw_plan_dft_r2c_3d(_grid_dimension[2], _grid_dimension[1], _grid_dimension[0], pot, pot_fft, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Divide by nr grid points 
    for (size_t i = 0 ; i<nr_grid_points_trans; i++){
        pot_fft[i][0] = pot_fft[i][0] / _nr_grid_points;
        pot_fft[i][1] = pot_fft[i][1] / _nr_grid_points;
    }

    //Construct incremental reciprocal vectors
    std::vector<std::vector<double>> reciprocal_vecs_incr;
    reciprocal_vecs_incr.push_back(cross_product(_lattice_vecs[1],_lattice_vecs[2],2.0*M_PI/_cell_volume));
    reciprocal_vecs_incr.push_back(cross_product(_lattice_vecs[2],_lattice_vecs[0],2.0*M_PI/_cell_volume));
    reciprocal_vecs_incr.push_back(cross_product(_lattice_vecs[0],_lattice_vecs[1],2.0*M_PI/_cell_volume));

    //Construct k-point mesh, consider factor due to redundancies in FT 
    std::vector<kpoint_t> kpoints = construct_kpoints(reciprocal_vecs_incr);

    ////////////////////////////////////
    // Do reciprocal space integration//
    ////////////////////////////////////

    // Do integration on reciprocal grid if basis functions on same centers and criterion for exponents fullfilled

    // Loop over shells
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nshell_prim,nr_grid_points_trans,kpoints,decontr_basis,decontr_shells,ints_prim,pot_fft,cout)
#endif

    for (size_t a = 0; a < nshell_prim ; a++){
        for (size_t b = 0; b <= a ; b++){
            if (decontr_shells[a].ind_center != decontr_shells[b].ind_center) {
                continue; // done in real space
            }
            else if (decontr_shells[a].exp[0]+decontr_shells[b].exp[0] <= _alpha_cut) {
                continue; // done in real space
            }
            else {
                // Initialize matrix for integrals over this shell pair
                arma::mat shell_ints(decontr_shells[a].ang_mom_comp.size(),decontr_shells[b].ang_mom_comp.size());
                shell_ints.zeros();

                int L = decontr_shells[a].ang_mom + decontr_shells[b].ang_mom;
                double alpha = decontr_shells[a].exp[0] + decontr_shells[b].exp[0];
                std::vector<tiger::ang_mom_comp_t> compA = decontr_shells[a].ang_mom_comp;
                std::vector<tiger::ang_mom_comp_t> compB = decontr_shells[b].ang_mom_comp;
                int multA = compA.size();
                int multB = compB.size();

                double factor = pow(sqrt(M_PI/alpha),3);
                double expfactor = -0.25 / alpha;

                // loop over kpoints
                size_t first = 0;
                for(size_t i=0; i<nr_grid_points_trans; i++, first++){
                    double x = kpoints[i].x;
                    double y = kpoints[i].y;
                    double z = kpoints[i].z;


                    // Determine phase factor exp(i*k*R)
                    double kR = x * decontr_shells[a].center.x  + y * decontr_shells[a].center.y  + z * decontr_shells[a].center.z ;
                    std::complex<double> phase(cos(kR),sin(kR));

                    // Determine integrals over gaussians 
                    std::vector<complex<double>> int_gauss_x = int_gauss_pw(x,L,alpha);
                    std::vector<complex<double>> int_gauss_y = int_gauss_pw(y,L,alpha);
                    std::vector<complex<double>> int_gauss_z = int_gauss_pw(z,L,alpha);

                    double prefactor = kpoints[i].factor * exp(expfactor *  (x*x + y*y + z*z)) ;

                    // calculate the final integrals
                    // loop over cartesian components of shell
                    for (size_t ai = 0; ai < multA ; ai++) {
                        for (size_t bi = 0; bi < multB ; bi++) {
                        	shell_ints(ai,bi) +=  prefactor * real(complex<double>(pot_fft[i][0],pot_fft[i][1]) * phase *
                                    int_gauss_x[compA[ai].l+compB[bi].l] *
                                    int_gauss_y[compA[ai].m+compB[bi].m] *
                                    int_gauss_z[compA[ai].n+compB[bi].n]);
                        }
                    }
                }

                for (size_t ai = 0; ai < multA ; ai++) {
                	for (size_t bi = 0; bi < multB ; bi++) {
                        ints_prim(decontr_shells[a].first_ind+ai , decontr_shells[b].first_ind+bi)=shell_ints(ai,bi)
                        		* decontr_shells[a].relnorms[ai] * decontr_shells[b].relnorms[bi]
                        		* decontr_shells[a].coeff_norm[0] * decontr_shells[b].coeff_norm[0]
                        		* factor ;
                    }
                }
            }
        }
    }

    // free fft memory
    fftw_free(pot_fft);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(nbas_prim,ints_prim)
#endif

    // add remaining elements based on symmetry: int(b,a)=int(a,b)
    for (size_t a=0; a < nbas_prim; a++) {
        for (size_t b =0; b < a; b++) {
            ints_prim(b,a)=ints_prim(a,b);
        }
    }


    return ints_prim;

}

void EmbeddingData::print_tigerci(arma::mat ints, std::string filename) {
    ofstream int_file;
    int_file.open(filename);
    int_file.precision(12);
    ints.raw_print(int_file);
    int_file.close();
}

void EmbeddingData::print_pyscf(arma::mat ints, std::string filename) {
    ofstream int_file;
    int_file.open(filename);
    int_file.precision(12);
    ints.raw_print(int_file);
    int_file.close();
}


void EmbeddingData::print_molpro(arma::mat ints, std::string filename) {
    ofstream int_file;
    int_file.open(filename);
    int_file.precision(12);
    int_file << "BEGIN_DATA," << endl << "# MATRIX EMB_INT            EMB    SYMMETRY=1" << endl;

    for (size_t i = 0; i!=ints.n_cols ; i++) {
        size_t n = 1;
        for (size_t j = 0; j!=ints.n_rows ; j++) {
            int_file << "    " << std::setprecision(12) << ints(j,i) << ",";
            if (n==5) {
                int_file << endl;
                n=1;
            }
            else n++;
        }
        if (n!=1){
        	int_file << endl;
        }
    }
    int_file << "END_DATA," << endl;
    int_file.close();
}

void EmbeddingData::print_gamess(arma::mat ints, std::string filename) {
    ofstream int_file;
    int_file.open(filename);
    int_file.precision(12);

    for (size_t i = 0; i!=ints.n_cols ; i++) {
        for (size_t j = 0; j<=i ; j++) {
            int_file << "    " << std::setprecision(12) << ints(j,i) << endl;
        }
    }
    int_file.close();
}

void EmbeddingData::print_molcas(arma::mat ints, std::string filename) {
    ofstream int_file;
    int_file.open(filename);

    BasisData *decontr_basis = _basis->decontract();
    std::vector<tiger::shell_t> shells = decontr_basis->get_shells();

    int natoms = shells.back().ind_center + 1;
    std::vector<int> atoms_ang_mom;
    // for every atom determine max angular momentum to be able to loop correctly when iterating over shells
	for (size_t i = 0; i<natoms; i++) {
		int ang_mom = 0;
		for (std::vector<tiger::shell_t>::iterator it = shells.begin(); it!=shells.end(); it++) {
			if (it->ind_center == i) {
				ang_mom = it->ang_mom;
			}
		}
		atoms_ang_mom.push_back(ang_mom);
	}

    // Molcas reader needs blocks for every pair of atoms separately -> iterate over atom pairs and search for corresponding basis functions
    for (size_t iA = 0; iA < natoms; iA++) {
    	for (size_t iB = 0; iB <= iA; iB++) {

    		size_t lB_max = atoms_ang_mom[iB]; // upper limit for running index over angular momenta on atom B; if A==B: lB <= lA
    		// Molcas wants to have for each atom pair one block for every pair of angular momenta l
    		for (size_t lA = 0; lA <= atoms_ang_mom[iA]; lA++) {
    			if (iA == iB) {
    				lB_max = lA;
    			}
    			for (size_t lB = 0; lB <= lB_max; lB++) {

    				// Molcas wants on the next level the full block of results for each m-component set
    				for (size_t mB = 0; mB < ((lB+1)*(lB+2)/2); mB++) {
    					for (size_t mA = 0; mA < ((lA+1)*(lA+2)/2); mA++) {
    						for (std::vector<tiger::shell_t>::iterator itB = shells.begin(); itB!=shells.end(); itB++) {
    							if (itB->ind_center != iB || itB->ang_mom!=lB) continue;
    							for (std::vector<tiger::shell_t>::iterator itA = shells.begin(); itA!=shells.end(); itA++) {
    								if (itA->ind_center != iA || itA->ang_mom!=lA) continue;
    								int_file.precision(3);
    								int_file << fixed << " embint2: la=" << setw(3) << itA->ang_mom << " lb=" << setw(3) << itB->ang_mom
    								    << " expA="<< setw(12) << itA->exp[0] << " expB" << setw(12) << itB->exp[0]
    								    << " Apos="<< setw(7) << itA->center.x << setw(7) << itA->center.y << setw(7) << itA->center.z
    								    << " RBpos="<< setw(7) << itB->center.x << setw(7) << itB->center.y << setw(7) << itB->center.z << " =>";
    								int_file.precision(12);
    								int_file << setw(24) << ints(itA->first_ind+mA,itB->first_ind+mB);
    								int_file << endl;
    							}
    						}
    					}
    				}
    			}
    		}
    	}
    }
    int_file.close();
}

arma::mat EmbeddingData::sort_molpro(arma::mat ints) {

    std::vector<tiger::shell_t> shells  = _basis->get_shells();
    std::vector<int> map;

    // iterate through shells and build the map according to the angular momenta
    for (std::vector<tiger::shell_t>::iterator it = shells.begin(); it!= shells.end(); it++) {
        int ang_mom = it->ang_mom;
        int offset = map.size();

        switch (ang_mom) {
            case 0: {
                for (size_t k = 0 ; k!=(ang_mom*2+1) ; k++)  map.push_back(_map_s[k] + offset);
                break;
            }
            case 1: {
                for (size_t k = 0 ; k!=(ang_mom*2+1) ; k++)  map.push_back(_map_p[k] + offset);
                break;
            }
            case 2: {
                for (size_t k = 0 ; k!=(ang_mom*2+1) ; k++)  map.push_back(_map_d[k] + offset);
                break;
            }
            case 3: {
                for (size_t k = 0 ; k!=(ang_mom*2+1) ; k++)  map.push_back(_map_f[k] + offset);
                break;
            }
            case 4: {
                for (size_t k = 0 ; k!=(ang_mom*2+1) ; k++)  map.push_back(_map_g[k] + offset);
                break;
            }
            case 5: {
                for (size_t k = 0 ; k!=(ang_mom*2+1) ; k++)  map.push_back(_map_h[k] + offset);
                break;
            }
            case 6: {
                for (size_t k = 0 ; k!=(ang_mom*2+1) ; k++)  map.push_back(_map_i[k] + offset);
                break;
            }
            default:  throw runtime_error("maximum angular momentum for embedding with Molpro is 6");
        }
    }

    arma::mat ints_sorted(ints.n_rows,ints.n_cols,arma::fill::zeros) ;

    for (size_t i = 0; i!=ints_sorted.n_cols ; i++) {
        for (size_t j=0; j!=ints_sorted.n_rows ; j++) {
            ints_sorted(map[j],map[i])=ints(j,i);
        }
    }

    return ints_sorted;
}

arma::mat EmbeddingData::sort_gamess(arma::mat ints) {

    std::vector<tiger::shell_t> shells  = _basis->get_shells();
    std::vector<int> map;

    // iterate through shells and build the map according to the angular momenta
    for (std::vector<tiger::shell_t>::iterator it = shells.begin(); it!= shells.end(); it++) {
        int ang_mom = it->ang_mom;
        int offset = map.size();
        int mult = (ang_mom+1)*(ang_mom+2)/2;

        switch (ang_mom) {
            case 0: {
                for (size_t k = 0 ; k!=mult ; k++)  map.push_back(_map_gamess_s[k] + offset);
                break;
            }
            case 1: {
                for (size_t k = 0 ; k!=mult ; k++)  map.push_back(_map_gamess_p[k] + offset);
                break;
            }
            case 2: {
                for (size_t k = 0 ; k!=mult ; k++)  map.push_back(_map_gamess_d[k] + offset);
                break;
            }
            case 3: {
            	for (size_t k = 0 ; k!=mult ; k++)  map.push_back(_map_gamess_f[k] + offset);
            	break;
            }
            default:  throw runtime_error("maximum angular momentum for embedding with GAMESS is 3");
        }
    }

    arma::mat ints_sorted(ints.n_rows,ints.n_cols,arma::fill::zeros) ;

    for (size_t i = 0; i!=ints_sorted.n_cols ; i++) {
        for (size_t j=0; j!=ints_sorted.n_rows ; j++) {
            ints_sorted(map[j],map[i])=ints(j,i);
        }
    }

    return ints_sorted; 
}

