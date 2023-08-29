/* 
 * File:   ConvertKeys.cpp
 * Author: Francis Ricci
 *
 * Created on June 10, 2014, 11:00 AM
 */

#include <string.h>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <memory>
#include <stdexcept>
#include <set>

#include "ConvertKeys.h"
#include "Keyword.h"
#include "InputReader.h"

// Unordered maps for the storage of name-value pairs.
// Must be global to allow for calling by "C" routines from Fortran.
static unordered_map<string, int> global_ints;
static unordered_map<string, double> global_doubles;
static unordered_map<string, bool> global_bools;
static unordered_map<string, string> global_strings;
static unordered_map<string, vector<vector<int>>> global_int_vecs;
static unordered_map<string, vector<vector<double>>> global_double_vecs;

using namespace std;

// convert the keys in key_map into variable name-value pairs,
// store in an unordered map
void ConvertKeys::keys_to_globals(InputReader& key_map) {
    bool fullyInt = false;
    vector<string> names = key_map.getNames();
    for (string name : names){
        if (name == "NUM THREADS")
            global_ints["numThreads"] = key_map.getIntKeyword(name);
        else if (name == "CARTESIAN")
            global_bools["spherical"] = !key_map.getBoolKeyword(name);
//TODO: Use this keyword
//        else if (name == "AO INTEGRAL THRESHOLD")
//            global_doubles["ao_integral_threshold"] = key_map.getDoubleKeyword(name);
        else if (name == "XYZ FILE"){
            my_strings["xyz_file"] = key_map.getStringKeyword(name);
        }
        else if (name == "BASIS FILE"){
            my_strings["basis_file"] = key_map.getStringKeyword(name);
        }
        else if (name == "GRID DIMENSIONS"){
            global_int_vecs["grid dimensions"] = key_map.getIntVecsKeyword(name);
        }
        else if (name == "LATTICE VECTORS"){
            global_double_vecs["lattice vectors"] = key_map.getDoubleVecsKeyword(name);
        }
        else if (name == "SHIFT VECTOR"){
            global_double_vecs["shift vector"] = key_map.getDoubleVecsKeyword(name);
        }
        else if (name == "EMBEDDING POTENTIAL FILE"){
            global_strings["embedding potential file"] = key_map.getStringKeyword(name);
        }

        else if (name == "EMBEDDING INTS REAL"){
            global_bools["embedding ints real"] = key_map.getBoolKeyword(name);
        }

        else if (name == "OUTPUT FORMAT"){
            global_strings["output format"] = key_map.getStringKeyword(name);
        }


        else
            throw runtime_error("Tried to read variable " + name +
                    " but it cannot be found");
    }
    // take care of a few reference things
    check_input();
}

double ConvertKeys::get_global_double(string name) const{
    return global_doubles.at(name);
}

int ConvertKeys::get_global_int(string name) const{
    return global_ints.at(name);
}

bool ConvertKeys::get_global_bool(string name) const{
    return global_bools.at(name);
}

string ConvertKeys::get_global_string(string name) const{
    return global_strings.at(name);
}

vector<vector<int>> ConvertKeys::get_global_int_vecs(string name) const{
    return global_int_vecs.at(name);
}

vector<vector<double>> ConvertKeys::get_global_double_vecs(string name) const{
    return global_double_vecs.at(name);
}

void ConvertKeys::add_global_int(string name, int val){
    global_ints[name] = val;
}

void ConvertKeys::add_global_bool(string name, bool val){
    global_bools[name]=val;
}

// get string val for global variable c_name
string ConvertKeys::get_my_string(string name) const{
    return my_strings.at(name);
}

// do some input parameter checking
void ConvertKeys::check_input(){

  // check number of threads
  int numThreads = global_ints.at("numThreads");
#if !defined(_OPENMP)
  if (numThreads > 1)
    throw runtime_error("ERROR: This is a serial code and you want more than one thread!");
#endif

  // check whether output format is available
  string method = global_strings.at("output format");
  std::set<std::string> methods = {"TIGERCI","MOLPRO","GAMESS","MOLCAS","NWCHEM","QCHEM","PYSCF"};
  if ( methods.find(method) == methods.end() ) {
    throw runtime_error("ERROR: This option for OUTPUT FORMAT is not implemented!");    
  }
}

