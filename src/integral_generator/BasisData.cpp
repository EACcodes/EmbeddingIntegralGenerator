/*
 * @file BasisData.cpp
 *
 *  Created on: Jan 15, 2015
 *      Author: Francis Ricci
 *
 *  Caroline Krauter:
 *  Added implementation of BasisData that is independent of Erkale
 */

#include "BasisData.h"
#include "EmbeddingData.h"

#include <map>

// pointer to an BasisData object (stored by main.cpp), which allows the
// Fortran interface functions to access the class objects (a workaround for
// Fortran not having access to C++ objects)
static BasisData *data;
static int64_t nCho;
static bool init = false;

/* C/Fortran interface functions */

void init_bas(BasisData *bas_data){
    if (init)
        return;
    data = bas_data;
    init = true;
}

BasisData::~BasisData(){}

EmbeddingBasis::EmbeddingBasis(ConvertKeys& global_vars) {
    _globals = &global_vars;
    _atoms = read_structure(_globals->get_my_string("xyz_file"));
    create_basis();
}


EmbeddingBasis::~EmbeddingBasis(){

}

BasisData* EmbeddingBasis::duplicate(){
	BasisData *new_bas = new EmbeddingBasis(*this, false);
    return new_bas;
}

BasisData* EmbeddingBasis::decontract(){
	BasisData *new_bas = new EmbeddingBasis(*this, true);
    return new_bas;
}

EmbeddingBasis::EmbeddingBasis(EmbeddingBasis& bas, bool decontracted){
    _globals = bas._globals;
    _atoms = bas._atoms;
    _trans = bas._trans; // assumes the contraction to be done wrt Cartesian functions

    if (!decontracted) {
        _shells = bas._shells;
    }

	else{ // decontract basis
		size_t index_shell = 0;
		int last_ang_mom  = -1;
		for (size_t i = 0; i<bas._shells.size(); i++){
			if (last_ang_mom == bas._shells[i].ang_mom){
				continue;
			}
			else {
				tiger::shell_t old_shell = bas._shells[i];
				last_ang_mom = old_shell.ang_mom;
				for (size_t j =0; j<old_shell.exp.size(); j++){
					tiger::shell_t new_shell;

					new_shell.center = old_shell.center;
					new_shell.ind_center = old_shell.ind_center;
					new_shell.ang_mom = old_shell.ang_mom;

					// always use Cartesian functions internally
					new_shell.use_sph  = false;
					new_shell.nbf = (new_shell.ang_mom +1)*(new_shell.ang_mom +2)/2;

					new_shell.coeff.push_back(1.0);
					// normalization coefficient
					if (_globals->get_global_string("output format")=="MOLCAS") {
						// Unnormalized functions in MOLCAS
						new_shell.coeff_norm.push_back(1.0);
						for (size_t k=0; k < old_shell.relnorms.size(); k++) {
							new_shell.relnorms.push_back(1.0);
						}
					}
					else {
						double factor = 1.0 / pow(2*old_shell.exp[j],old_shell.ang_mom+1.5) * pow(M_PI,1.5)
									/ pow(2,old_shell.ang_mom) * dfact(2*old_shell.ang_mom-1) ;
						new_shell.coeff_norm.push_back(1.0 / sqrt(factor));
						new_shell.relnorms  = old_shell.relnorms;
					}

					new_shell.exp.push_back(old_shell.exp[j]);
					new_shell.first_ind = index_shell;

					new_shell.ang_mom_comp = old_shell.ang_mom_comp;

					index_shell+=new_shell.nbf;
					_shells.push_back(new_shell);
				}
			}
		}
    }
   //print_basis();
}


int EmbeddingBasis::get_Ncart(){
	int nBas = 0;
	for (size_t i=0; i<_shells.size(); i++){
		nBas += (_shells[i].ang_mom + 2)*(_shells[i].ang_mom + 1)/2;
	}
	return nBas;
}

int EmbeddingBasis::get_Nsph(){
	int nBas = 0;
	for (size_t i=0; i<_shells.size(); i++){
		nBas += 2 * _shells[i].ang_mom + 1;
	}
	return nBas;
}


std::vector<tiger::coords_t> EmbeddingBasis::get_coordinates(){

	std::vector<tiger::coords_t> coordinates;

	for (std::vector<tiger::atom_t>::iterator it = _atoms.begin(); it< _atoms.end(); it++){
		coordinates.push_back(it->coords);
	}

	return coordinates;
}

std::vector<tiger::atom_t> EmbeddingBasis::get_atoms(){
	return _atoms;
}

std::vector<int> EmbeddingBasis::get_nFuncInShell(){
	std::vector<int> funcInShell;
	for (size_t ind = 0; ind < _shells.size(); ind++) {
			funcInShell.push_back(_shells[ind].nbf);
	}
	if (funcInShell.size()!= get_nShell()){
		throw runtime_error("Error in counting functions in Shells.");
	}
	return funcInShell;
}

int EmbeddingBasis::get_nBas(){
	std::vector<int> funcInShell = get_nFuncInShell();
	int nBas = 0;
	for (size_t ind = 0; ind < funcInShell.size(); ind++) {
		nBas+=funcInShell[ind];
	}
	return nBas;
}

int EmbeddingBasis::get_nShell(){
	return _shells.size();
}

int EmbeddingBasis::get_nAtoms(){
	return _atoms.size();
}

std::vector<int> EmbeddingBasis::get_basis_map(){
	// currently not used
}

arma::vec EmbeddingBasis::eval_func(double x, double y, double z){

	int nbas_prim = get_nBas();
    arma::vec ints(nbas_prim);
    ints.zeros();

    int index_vec = 0;
    for (int ind =0 ; ind < get_nShell() ; ind++) {
    	arma::vec shell_ints = eval_func(ind,x,y,z);

    	// insert in vector of ints
         for (arma::vec::iterator it2 = shell_ints.begin(); it2 != shell_ints.end() ; it2++,index_vec++) {
        	 ints(index_vec) = *it2;
         }
    }
    return ints;
}

arma::vec EmbeddingBasis::eval_func(size_t shell_ind, double x, double y, double z){

	double alpha = _shells[shell_ind].exp[0];
	double coeff_norm = _shells[shell_ind].coeff_norm[0];
	tiger::coords_t center = _shells[shell_ind].center;
	double exp1 = alpha * ((x - center.x)*(x - center.x) + (y - center.y)*(y - center.y) + (z - center.z)*(z - center.z));
	std::vector<tiger::ang_mom_comp_t> ang_mom_comp = _shells[shell_ind].ang_mom_comp;

	int i=0;
	arma::vec shell_ints(ang_mom_comp.size());
	shell_ints.zeros();
	for (std::vector<tiger::ang_mom_comp_t>::iterator it = ang_mom_comp.begin() ; it != ang_mom_comp.end() ; it++, i++) {
		shell_ints(i) = coeff_norm * _shells[shell_ind].relnorms[i] * pow(x-center.x,it->l) * pow(y-center.y,it->m) * pow(z-center.z,it->n) * exp(-exp1);
	}

    return shell_ints;
}

void EmbeddingBasis::print_basis() {
	cout << endl;
	cout << "#######################################" << endl;
	cout << "#######  Basis set information  #######" << endl;
	cout << "#######################################" << endl << endl;

	for (size_t i=0; i<_shells.size();i++){
		print_shell(i);
	}

	cout << "#######################################" << endl;
}

void EmbeddingBasis::print_shell(size_t ind) {
	cout << endl << "Shell with index " << ind << " at position " << _shells[ind].center.x
			<< " " << _shells[ind].center.y << " " << _shells[ind].center.z <<  endl;
	cout << "   center   ang. mom.   ang. mom. comp.      exp.      coeff.   norm. coeff.   rel. norm coeff." << endl;

	for (size_t i=0; i<_shells[ind].ang_mom_comp.size(); i++){
		cout << "      " << setw(11) << _shells[ind].ind_center << setw(11) <<  _shells[ind].ang_mom
				<<  "(" << _shells[ind].ang_mom_comp[i].l << ","
				<< _shells[ind].ang_mom_comp[i].m << "," << _shells[ind].ang_mom_comp[i].n << ")";
		for (size_t n=0; n < _shells[ind].exp.size(); n++) {
			if (n==0){
				cout.precision(4);
				cout << fixed << "       " << setw(12) << _shells[ind].exp[n]
				     << setw(12) << _shells[ind].coeff[n] << setw(16) << _shells[ind].coeff_norm[n] << _shells[ind].relnorms[i] <<  endl;
			}
			else {
				cout << "                                          " << setw(12) <<  _shells[ind].exp[n]
				     << setw(12) << _shells[ind].coeff[n] << setw(16) << _shells[ind].coeff_norm[n] <<  endl;
			}
		}
	}
}

std::vector<tiger::shell_t> EmbeddingBasis::get_shells() {
	return _shells;
}

arma::mat EmbeddingBasis::get_trans_mat() {
	return _trans;
}

bool EmbeddingBasis::is_spherical(size_t ind) {
	return _shells[ind].use_sph;
}

arma::mat EmbeddingBasis::get_trans_cart_sph(size_t ind) {
	return trans_cart_sph(_shells[ind].ang_mom);
}

arma::mat EmbeddingBasis::get_trans_cart_sph() {
	arma::mat trans(get_Nsph(),get_Ncart());
	trans.zeros();

	size_t ind_sph = 0;
	size_t ind_cart = 0;

	for (size_t i=0; i < _shells.size(); i++){
		int l = _shells[i].ang_mom;
		arma::mat trans_shell;
		if (l>1){
			trans_shell = trans_cart_sph(l);
		}
		else {
			trans_shell.resize(2*l+1, 2*l+1);
			trans_shell.zeros(2*l+1, 2*l+1);
			for (size_t m=0; m<(2*l+1); m++){
				trans_shell(m,m) = 1.0;
			}
		}
		for (size_t j=0; j<trans_shell.n_rows; j++){
			for (size_t k=0; k<trans_shell.n_cols; k++){
				trans(ind_sph+j,ind_cart+k) = trans_shell(j,k);
			}
		}
		ind_cart += _shells[i].ang_mom_comp.size();
		ind_sph  += 2*_shells[i].ang_mom +1;
	}
	// sanity check for dimensions
	if (ind_cart != get_Ncart() || ind_sph != get_Nsph()){
		cout << "Ncart: " << get_Ncart() << " and Nsph: " << get_Nsph() << " but indices: " << ind_cart << " " << ind_sph << endl;
		throw runtime_error("Clash of dimensions in forming transformation matrix between Cartesian and spherical coordinates");
	}
	return trans;
}

size_t EmbeddingBasis::get_shell_center_ind(size_t ind) {
	return _shells[ind].ind_center;
}

void EmbeddingBasis::create_basis(){

	// list to hold lines after deleting comments, converting to lower case and deleting empty lines
    std::list<string> lines_basis = read_basis_file(_globals->get_my_string("basis_file"));

    size_t ind_counter = 0; // counter to keep track of tiger-type shells

    //look for basis set for each atom and convert to internal format
    for (std::vector<tiger::atom_t>::iterator it = _atoms.begin(); it < _atoms.end(); it++) {
    	// Read molcas format
    	std::vector<tiger::molcas_shell_t> shells = read_shells(it->Z, lines_basis);

    	// Convert to internal definition of shells
    	for (std::vector<tiger::molcas_shell_t>::iterator it2 = shells.begin(); it2!=shells.end(); it2++){
    		for (size_t i=0; i < it2->n_contr; i++){
    			tiger::shell_t shell;
    			shell.center = it->coords;
    			shell.ind_center = std::distance(_atoms.begin(), it);
    			shell.ang_mom = it2->ang_mom;
    			// always use Cartesian basis functions internally
				shell.nbf = (shell.ang_mom +1)*(shell.ang_mom +2)/2;
				shell.use_sph  = false;

    			// setup Cartesian angular momentum components
    			for (size_t a = 0; a <= shell.ang_mom ; a++){
    				for (size_t b = 0; b <= a ; b++) {
    					tiger::ang_mom_comp_t ang_mom;
    					ang_mom.l = shell.ang_mom - a;
    					ang_mom.m = a - b;
    					ang_mom.n = b;
    					shell.ang_mom_comp.push_back(ang_mom);
    				}
    			}
    			// sanity check for length of ang_mom_comp
    			if (shell.ang_mom_comp.size() != shell.nbf) {
    				throw runtime_error("Error in basis functions setup. Unexpected number of angular momentum components in a shell.");
    			}

    			// setup exponents and coefficients
    			// TODO: Only use primitive function if coefficient 0
    			shell.exp = it2->exp;
    			for (size_t n = 0; n < shell.exp.size(); n++){
    				shell.coeff.push_back(it2->contr_mat(n,i));
    			}
    			// normalize coefficients. Part which is not angular momentum component dependent
    			double factor = 0.0;
    			for (size_t ni = 0; ni < it2->exp.size(); ni++){
    				for (size_t nj = 0; nj < it2->exp.size(); nj++){
    					factor += shell.coeff[ni]*shell.coeff[nj]/pow(shell.exp[ni]+shell.exp[nj],shell.ang_mom+1.5);
    			    }
    			}
    			factor *= pow(M_PI,1.5) / pow(2,shell.ang_mom);
    			factor = 1.0 / sqrt(factor);
    			for (size_t n = 0; n < shell.coeff.size(); n++){
					shell.coeff_norm.push_back(factor*shell.coeff[n]);
    			}

    			// determine relative norms of different ang mom components
    			for (size_t n = 0; n < shell.ang_mom_comp.size(); n++){
    				shell.relnorms.push_back(sqrt( dfact(2*shell.ang_mom -1 ) / ( dfact(2*shell.ang_mom_comp[n].l - 1)
    						* dfact(2*shell.ang_mom_comp[n].m - 1) * dfact(2*shell.ang_mom_comp[n].n - 1))));
    			}

    			shell.first_ind = ind_counter;
    			ind_counter += shell.nbf;

    			_shells.push_back(shell);


    		}
    		// determine contraction matrix
   			int n_rows = _trans.n_rows;
   			int n_cols = _trans.n_cols;

   			//double mult = (it2->ang_mom +1)*(it2->ang_mom +2)/2;
   			double mult = (it2->ang_mom +1)*(it2->ang_mom +2)/2;;
   			_trans.resize(n_rows + (it2->n_prim * mult), n_cols + (it2->n_contr * mult));

   			for (size_t i=0; i < it2->n_prim; i++){
   				for (size_t j=0; j < it2->n_contr; j++){
   					for (size_t m=0; m < mult; m++){
   						_trans(n_rows+i*mult+m,n_cols+j*mult+m) = it2->contr_mat(i,j);
   					}
   				}
   			}
    	}
    }
    //print_basis();
}

std::list<string> EmbeddingBasis::read_basis_file(string file_name) {

	ifstream inputf;
	inputf.open(file_name);
	if (!inputf){
		throw runtime_error("Problem opening basis set file.");
	}

	std::list<string> lines;

	string line;
	while (getline(inputf, line)) {
		size_t found_comment = line.find('!');

		if (found_comment != std::string::npos) {
			// comments begin with ! or *
			if (line[0]=='!') {
				// whole line is comment -> skip
				continue;
			}
			else {
				// truncate line to delete comment
				line = line.substr(0,found_comment);
			}
		}

		found_comment = line.find('*');
		if (found_comment != std::string::npos) {
			// comments begin with ! or *
			if (line[0]=='!') {
				// whole line is comment -> skip
				continue;
			}
			else {
				// truncate line to delete comment
				line = line.substr(0,found_comment);
			}
		}
		// convert to lower case
		char c;
		int i=0;

		//cout << "string original: " << line << endl;

		while(line[i]){
			c=line[i];
			line[i]=tolower(c);
			i++;
		}

		//cout << "string lower: " << line << endl;

		// check whether line contains anything and delete empty lines
		std::stringstream lineStream(line);
		std::string token;
		size_t n = 0;
		while(lineStream >> token) {
			n++;
		}
		if (n==0){
			continue;
		}
		else {
			lines.push_back(line);
		}
	}
	return lines;
}

std::vector<tiger::molcas_shell_t> EmbeddingBasis::read_shells(double nuc_charge, std::list<string> &lines){

	std::vector<tiger::molcas_shell_t> shells;

    int max_l = 0;
    for (std::list<string>::iterator it = lines.begin(); it!=lines.end(); it++){
    	size_t found_basis = it->find("basis set");
    	if (found_basis != std::string::npos){
    		// found a new basis set, see whether it is the correct one

    		double atomic_charge = 0.;
    		std::list<string>::iterator it2 = it;
    		std::advance(it2,2); // information about nuclear charge and maximum l is contained two rows below

    		std::stringstream lineStream(*it2);
    		std::string token;
    		size_t n = 0; // counter to count elements in this line
    		while(lineStream >> token) {
    			if (n==0){
    				atomic_charge = std::stod(token);
    				n++;
    			}
    			else if (n==1){
    				max_l = std::stoi(token);
    				n++;
    			}
    			else {
    				cout << "WARNING! in basis set file a line supposed to contain the atomic charge and maximum l "
    						"is containing more values than expected. Check the results." << endl;
    			}
    		}

    		if (atomic_charge == nuc_charge){
    			// read next line containing the number of primitive and contracted functions for l=0
    			for (size_t l = 0; l <= max_l; l++) {
    				tiger::molcas_shell_t shell;

    				shell.ang_mom = l;

    				it2++;
    				n=0;
    				std::stringstream lineStream(*it2);
    				while(lineStream >> token) {
    					if (n==0){
    						shell.n_prim = std::stoi(token);
    						n++;
    					}
    					else if (n==1){
    						shell.n_contr = std::stoi(token);
    						n++;
    					}
    					else {
    						cout << "WARNING! in basis set file a line supposed to contain the number of "
    								"primitive and contracted functions for atom with charge "
    								<< nuc_charge << " and l=" << l <<
    								" contained more elements than expected. Check the results. " << endl;
    					}
    				}
    				// read exponents for this l
    				for (size_t e =0; e < shell.n_prim; e++){
    					it2++;
    					n=0;
    					std::stringstream lineStream(*it2);
    					while(lineStream >> token) {
    						if (n==0){
    							shell.exp.push_back(std::stod(token));
    						}
    						else {
    							cout << "WARNING! in basis set file a line supposed to contain only an exponent for atom with charge "
    						    	<< nuc_charge << " and l=" << l << " contained more elements than expected. Check the results. " << endl;
    						}
    					}
    				}

    				// read contraction coefficients for this l
    				shell.contr_mat.zeros(shell.n_prim,shell.n_contr);
    				for (size_t e = 0; e < shell.n_prim; e++){
    					it2++;
    					n=0;
    					std::stringstream lineStream(*it2);
    					while(lineStream >> token) {
    						if (n < shell.n_contr){
    							shell.contr_mat(e,n) = std::stod(token);
    							n++;
    						}
    						else {
    							cout << "WARNING! in basis set file a line supposed to contain "
    									"only contraction coeffiecients for atom with charge "
    						    		<< nuc_charge << " and l=" << l <<
    						    		" contained more elements than expected. Check the results. " << endl;
    						}
    					}
    				}
    				// calculate overlap S and scale coefficients by sqrt(S)

    				for (size_t i = 0; i< shell.n_contr; i++) {
    					double S = 0 ;
    					for (size_t j = 0; j< shell.n_prim; j++) {
    						for (size_t k=0; k< shell.n_prim; k++) {
    							double sum_exp = shell.exp[j]+shell.exp[k];
    							S+= shell.contr_mat(j,i)*shell.contr_mat(k,i)*pow(4*shell.exp[j]*shell.exp[k]/sum_exp/sum_exp,(2.*shell.ang_mom+3.)/4.);
    						}
    					}
    					for (size_t j=0; j< shell.n_prim; j++) {
    						shell.contr_mat(j,i)*=1./sqrt(S);
    					}

    				}

    				shells.push_back(shell);
    			}
    		}
    	}
    }

    if (shells.size()==0) {
    	ostringstream oss;
    	oss << "Could not find basis set for atom with atomic charge " << nuc_charge ;
    	throw runtime_error(oss.str());
    }

    return shells;
}


std::vector<tiger::atom_t> EmbeddingBasis::read_structure(string file_name){

	std::vector<tiger::atom_t> atoms;
	size_t nr_atoms;
	bool bohr = false;

	ifstream xyzf;
	xyzf.open(file_name);
	if (!xyzf){
		throw runtime_error("Problem opening the structure file.");
	}

	// First line containes number of atoms
	std::string line;
	getline(xyzf,line);
	std::stringstream lineStream(line);

	std::string token;
	size_t n = 0;
	while(lineStream >> token) {
		if (n==0){
			nr_atoms = std::stoi(token);
	    	n++;
	    }
	    else {
	    	cout << "WARNING! in xyz set file the first line is supposed to contain only the number of atoms."
	    			"Detected more elements than expected." << endl;
	    }
	}

	// Second line is a comment, unless "bohr" is contained (-> adjust unit; default: angstrom)
	getline(xyzf,line);
	// Convert line to lower case
	int i=0;
	char c;
	while(line[i]){
		c=line[i];
		line[i]=tolower(c);
		i++;
	}
	size_t found_bohr = line.find("bohr");
	if (found_bohr != std::string::npos) {
		bohr = true;
	}

	tiger::atom_t atom;
	// Third to nr_atoms+2 lines contain geometry data (default: angstrom)
	for (size_t i = 0; i<nr_atoms; i++){
		getline(xyzf,line);
		n = 0;
		std::stringstream lineStream(line);
		while(lineStream >> token) {
			if (n==0){
				atom.Z = nuc_chg(token);
		    	n++;
		    }
			else if (n==1) {
				atom.coords.x = std::stod(token);
				n++;
		    }
			else if (n==2) {
				atom.coords.y = std::stod(token);
				n++;
 			}
			else if (n==3) {
				atom.coords.z = std::stod(token);
				n++;
			}
			else {
				cout << "WARNING: a line in the xyz file specifying the geometry contains more elements than expected" << endl;
			}
		}

		if (!bohr){
			atom.coords.x = atom.coords.x * _ang2bohr;
			atom.coords.y = atom.coords.y * _ang2bohr;
			atom.coords.z = atom.coords.z * _ang2bohr;
		}
		atoms.push_back(atom);
	}
	return atoms;
}

int EmbeddingBasis::dfact(int number){

	int prod = 1.0;
	for (int t = number; t>=1; t-=2){
		prod*=t;
	}
	return prod;
}

int EmbeddingBasis::fact(int number){
	if (number < 0) {
		throw runtime_error("Trying to calculate factorial of a negative number!");
	}

	int factorial = 1;

	for (size_t i=1; i<=number; i++){
		factorial *= i;
	}

	return factorial;
}

int EmbeddingBasis::bin_coeff(int n, int k){
	int coeff;
	if (k<0 || k>n){
		coeff = 0;
	}
	else if (n == k) {
		coeff = 1;
	}
	else {
		coeff = fact(n)/(fact(k)*fact(n-k));
	}
	return coeff;
}

arma::mat EmbeddingBasis::trans_cart_sph(int l) {
	// according to Equ. 15 in Int.J.Quant.Chem. 54, 83-87, 1995.
	//    needs to be modified to real spherical harmonics to be useful
	//	  (Y^l_m + Y^l_-m)/sqrt(2) and (Y^l_m - Y^l_-m)/sqrt(-2)
	arma::mat trans(2*l +1, (l+1)*(l+2)/2);
	trans.zeros();

	for (size_t m = 0; m <= l; m++) {
		// iterate over Cartesian angular momentum components
		size_t ind_cart = 0;
		for (size_t a = 0; a <= l ; a++){
			for (size_t b = 0; b <= a ; b++,ind_cart++) {
				int lx = l - a;
				int ly = a - b;
				int lz = b;
				int j=(lx+ly-m)/2;

				double prod=0.;
				if ( (lx+ly-m)%2 == 0 && j>=0 ) {
					for (size_t i = 0; i <= (l-m)/2; i++) {
						prod += bin_coeff(l,i)*bin_coeff(i,j)*pow(-1,i)*fact(2*l-2*i)/(double)fact(l-m-2*i);
					}
					double factor = sqrt((double)fact(2*lx)*fact(2*ly)*fact(2*lz)*fact(l)*fact(l-m)
									/fact(2*l)/fact(lx)/fact(ly)/fact(lz)/fact(l+m))
									/pow(2,l)/fact(l);
					//cout << "Factor " << factor << " prod " << prod << endl;
					prod *= factor;


					//cout << "Result for lx=" << lx << ", ly=" << ly << ", lz=" << lz << " and m=" << m << " is " << prod << endl;

					int sum = 0;
					//cout << " j " << j << endl;
					for (size_t k=0; k<=j; k++) {
						if ( ((m-lx+2*k)/2)%2 == 0 ) {
							sum+=bin_coeff(j,k)*bin_coeff(m,lx-2*k);
						}
						else {
							sum-=bin_coeff(j,k)*bin_coeff(m,lx-2*k);
						}
					}


					if (m==0) {
						// m=0 is always real and same for complex and real harmonics
						trans(l,ind_cart) = prod * sum;
					}
					else if ((m-lx)%2 == 0) {
						// real component entering +m component of real spherical harmonics
						trans(m+l,ind_cart) = prod * sum * sqrt(2.);
					}
					else {
						// complex component entering -m component of real spherical harmonics
						trans(-m+l,ind_cart) = prod * sum * sqrt(2.);
					}

				}
			}
		}
	}
	return trans;
}

arma::mat EmbeddingBasis::trans_cart_sph_table(int l) {
	// according to Table1 in Int.J.Quant.Chem. 54, 83-87, 1995.
	// Maximum angular momentum available: 5
	// Equation 15 gives general equation for transformation coefficients
	//    need to be modified to real spherical harmonics to be useful

	arma::mat trans(2*l +1, (l+1)*(l+2)/2);
	trans.zeros();

	switch (l) {
		case 0: {
			trans(0,0)=1.;
			break;
		}
		case 1: {
			trans(0,1)=trans(1,2)=trans(2,0)=1.;
			break;
		}
		case 2: {
			trans(0,1)=trans(1,4)=trans(2,5)=trans(3,2)=1.;
			trans(2,0)=trans(2,3)=-0.5;
			trans(4,0)=sqrt(3)/2.;
			trans(4,3)=-sqrt(3)/2.;
			break;
		}
		case 3: {
			trans(0,1)=3.*sqrt(2.)/4.;
			trans(0,6)=-sqrt(10.)/4.;
			trans(1,4)=trans(3,9)=1.;
			trans(2,1)=trans(4,3)=-sqrt(6./5.)/4.;
			trans(2,6)=trans(4,0)=-sqrt(6.)/4;
			trans(2,8)=trans(4,5)= sqrt(6./5.);
			trans(3,2)=trans(3,7)=-3./2./sqrt(5.);
			trans(5,2)= sqrt(3.)/2.;
			trans(5,7)=-trans(5,2);
			trans(6,0)=-trans(0,6);
			trans(6,3)=-trans(0,1);
			break;
		}
		case 4: {
			trans(0,1)=sqrt(5.)/2.;
			trans(0,6)=-trans(0,1);
			trans(1,4)=3.*sqrt(2.)/4.;
			trans(1,11)=-sqrt(10.)/4.;
			trans(2,1)=trans(2,6)=-sqrt(5.)/2./sqrt(7.);
			trans(2,8)=3./sqrt(7.);
			trans(3,4)=trans(5,7)=-3.*sqrt(2.)/sqrt(7.)/4.;
			trans(3,11)=trans(5,2)=-3*sqrt(10.)/4./sqrt(7.);
			trans(3,13)=trans(5,9)=sqrt(10./7.);
			trans(4,0)=trans(4,10)=3./8.;
			trans(4,3)=3.*sqrt(3.)/4./sqrt(35.);
			trans(4,5)=trans(4,12)=-3.*sqrt(3.)/sqrt(35.);
			trans(4,14)=1.;
			trans(6,0)=-sqrt(5.)/4.;
			trans(6,5)=3.*sqrt(3.)/2./sqrt(7.);
			trans(6,10)=-trans(6,0);
			trans(6,12)=-trans(6,5);
			trans(7,2)=-trans(1,11);
			trans(7,7)=-trans(1,4);
			trans(8,0)=trans(8,10)=sqrt(35.)/8.;
			trans(8,3)=-3*sqrt(3.)/4.;
			break;
		}
		default: cout << "Trying to handle function with angular momentum "
				<< l << endl; throw runtime_error("This angular momentum "
						"not implemented in transformation from Cartesian "
						"to spherical functions.");
	}

	cout << "Transformation matrix Cart/Sph " << endl << trans << endl;
	return trans;
}


arma::mat EmbeddingBasis::trans_cart_sph_old(int l) {
	// according to https://en.wikipedia.org/wiki/Solid_harmonics

	// Originally used, but requires weird normalization of Cartesian basis functions
	// Could not reproduce reference values for l>3
	// switched to formulation according to Equ. 15 in Int.J.Quant.Chem. 54, 83-87, 1995.
	// left this function in case needed in the future
	arma::mat trans(2*l +1, (l+1)*(l+2)/2);
	trans.zeros();

	// due to Racah's normalization
	double prefact = sqrt( (2*l+1)/ (4*M_PI) );

	// m-independent prefactor from \Pi_l^m
	prefact *= pow (2.0,-l);

	for (int m=0; m<=l; m++){

		// prefactors in C_l^m and S_l^m; depend only on l and m
		double prefact_CS;
		if (m!=0){
			prefact_CS = sqrt(2.*fact(l-m)/fact(l+m)) * prefact;
		}
		else{
			prefact_CS = sqrt(1.*fact(l-m)/fact(l+m)) * prefact;
		}

		for (size_t k=0; k<=(l-m)/2.; k++){
			// prefactors of \Pi_l^m; m-independent part has been included above in prefact
			double prefact_pi = pow(-1,k) * bin_coeff(l,k) * bin_coeff(2*l-2*k,l) * prefact_CS;
			if (m!=0){
				prefact_pi*= fact(l-2*k)/fact(l-2*k-m);
			}

			for (size_t p=0; p<=m; p++){
				// common prefactor of A_m and B_m
				double prefact_AB = bin_coeff(m,p) * prefact_pi;

				// determine exponents of x,y,z
				for (size_t a=0; a<=k; a++){
					for (size_t b=0; b<=(k-a); b++){
						double factor = fact(k)/fact(a)/fact(b)/fact(k-a-b) * prefact_AB;
						int exp_x = p+2*a;
						int exp_y = m-p+2*b;
						int exp_z = l-m-2*a-2*b;
						int ind_cart = (exp_y+exp_z)*(exp_y+exp_z+1)/2 + exp_z;

						// \Pi_l^m * A_m and \Pi_l^m * B_m
						if (m > 0){
							int j = (m-p)%4;
							double A = 0.;
							double B = 0.;

							// determine cases when A or B != 0
							if (j == 0){ A = 1.;} // cos(0)
							else if (j == 1){ B = 1.;} // sin(\pi/2)
							else if (j == 2){ A = -1.;} // cos(\pi)
							else if (j == 3){ B = -1.;} // sind(3\pi/2)
							else {throw runtime_error("Error in calculating A_m(x,y) and B_m(x,y)");}

							trans( m+l,ind_cart) += factor * A;
							trans(-m+l,ind_cart) += factor * B;

						}
						else {
							trans(l,ind_cart) += factor;
						}
					}
				}
			}
		}
	}
	return trans;
}



int EmbeddingBasis::nuc_chg(std::string token){

	std::map<string,int> nuc_chg;

	nuc_chg["H" ] =    1;
	nuc_chg["HE"] =    2;
	nuc_chg["LI"] =    3;
	nuc_chg["BE"] =    4;
	nuc_chg["B" ] =    5;
	nuc_chg["C" ] =    6;
	nuc_chg["N" ] =    7;
	nuc_chg["O" ] =    8;
	nuc_chg["F" ] =    9;
	nuc_chg["NE"] =   10;
	nuc_chg["NA"] =   11;
	nuc_chg["MG"] =   12;
	nuc_chg["AL"] =   13;
	nuc_chg["SI"] =   14;
	nuc_chg["P" ] =   15;
	nuc_chg["S" ] =   16;
	nuc_chg["CL"] =   17;
	nuc_chg["AR"] =   18;
	nuc_chg["K" ] =   19;
	nuc_chg["CA"] =   20;
	nuc_chg["SC"] =   21;
	nuc_chg["TI"] =   22;
	nuc_chg["V" ] =   23;
	nuc_chg["CR"] =   24;
	nuc_chg["MN"] =   25;
	nuc_chg["FE"] =   26;
	nuc_chg["CO"] =   27;
	nuc_chg["NI"] =   28;
	nuc_chg["CU"] =   29;
	nuc_chg["ZN"] =   30;
	nuc_chg["GA"] =   31;
	nuc_chg["GE"] =   32;
	nuc_chg["AS"] =   33;
	nuc_chg["SE"] =   34;
	nuc_chg["BR"] =   35;
	nuc_chg["KR"] =   36;
	nuc_chg["RB"] =   37;
	nuc_chg["SR"] =   38;
	nuc_chg["Y" ] =   39;
	nuc_chg["ZR"] =   40;
	nuc_chg["NB"] =   41;
	nuc_chg["MO"] =   42;
	nuc_chg["TC"] =   43;
	nuc_chg["RU"] =   44;
	nuc_chg["RH"] =   45;
	nuc_chg["PD"] =   46;
	nuc_chg["AG"] =   47;
	nuc_chg["CD"] =   48;
	nuc_chg["IN"] =   49;
	nuc_chg["SN"] =   50;
	nuc_chg["SB"] =   51;
	nuc_chg["TE"] =   52;
	nuc_chg["I" ] =   53;
	nuc_chg["XE"] =   54;
	nuc_chg["CS"] =   55;
	nuc_chg["BA"] =   56;
	nuc_chg["LA"] =   57;
	nuc_chg["CE"] =   58;
	nuc_chg["PR"] =   59;
	nuc_chg["ND"] =   60;
	nuc_chg["PM"] =   61;
	nuc_chg["SM"] =   62;
	nuc_chg["EU"] =   63;
	nuc_chg["GD"] =   64;
	nuc_chg["TB"] =   65;
	nuc_chg["DY"] =   66;
	nuc_chg["HO"] =   67;
	nuc_chg["ER"] =   68;
	nuc_chg["TM"] =   69;
	nuc_chg["YB"] =   70;
	nuc_chg["LU"] =   71;
	nuc_chg["HF"] =   72;
	nuc_chg["TA"] =   73;
	nuc_chg["W" ] =   74;
	nuc_chg["RE"] =   75;
	nuc_chg["OS"] =   76;
	nuc_chg["IR"] =   77;
	nuc_chg["PT"] =   78;
	nuc_chg["AU"] =   79;
	nuc_chg["HG"] =   80;
	nuc_chg["TL"] =   81;
	nuc_chg["PB"] =   82;
	nuc_chg["BI"] =   83;
	nuc_chg["PO"] =   84;
	nuc_chg["AT"] =   85;
	nuc_chg["RN"] =   86;
	nuc_chg["FR"] =   87;
	nuc_chg["RA"] =   88;

	for (size_t i = 0; i< token.length() ; i++) {
		token[i] = toupper(token[i]);
	}
	return nuc_chg[token];
}































