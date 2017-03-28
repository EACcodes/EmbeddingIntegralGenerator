/*
 * File:   InputReader.cpp
 * Authors: Francis Ricci and David Krisiloff
 *
 * \todo Add sanity checks on different inputs
 */

#include <string>
#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include "Keyword.h"
#include "StringUtils.h"
#include "InputReader.h"

using namespace std;

InputReader::InputReader() {
    /*
     * This is where we define each of the input_keys in the TigerCI
     * input file, their default values, nonlocal values, etc..
     */

    // Tolerance and thresholds
//    input_keys["AO INTEGRAL THRESHOLD"] = shared_ptr<Keyword>
//    (new DoubleKeyword("AO INTEGRAL THRSHOLD", 1.0e-12));

    // Flags
    input_keys["CARTESIAN"] = shared_ptr<Keyword>
    (new FlagBoolKeyword("CARTESIAN"));
    input_keys["EMBEDDING INTS REAL"] = shared_ptr<Keyword>
    (new FlagBoolKeyword("EMBEDDING INTS REAL"));

    // Other parameters
    input_keys["NUM THREADS"] = shared_ptr<Keyword>
    (new DefaultIntKeyword("NUM THREADS", 1));

    //Vector parameters
    input_keys["GRID DIMENSIONS"] = shared_ptr<Keyword>
    (new RequiredIntVecsKeyword("GRID DIMENSIONS"));
    input_keys["LATTICE VECTORS"] = shared_ptr<Keyword>
    (new RequiredDoubleVecsKeyword("LATTICE VECTORS"));
    input_keys["SHIFT VECTOR"] = shared_ptr<Keyword>
    (new DefaultDoubleVecsKeyword("SHIFT VECTOR",{{0.0,0.0,0.0}}));

    // Files and directories

    // Stand-alone related variables
    input_keys["XYZ FILE"] = shared_ptr<Keyword>
    (new RequiredStringKeyword("XYZ FILE"));
    input_keys["BASIS FILE"] = shared_ptr<Keyword>
    (new RequiredStringKeyword("BASIS FILE"));
    input_keys["EMBEDDING POTENTIAL FILE"] = shared_ptr<Keyword>
    (new DefaultStringKeyword("EMBEDDING POTENTIAL FILE","embpot.dat"));
    input_keys["OUTPUT FORMAT"] = shared_ptr<Keyword>
    (new RequiredStringKeyword("OUTPUT FORMAT"));
}

// Read and store input parameters from string file_name
void InputReader::read_tiguar_input(string file_name) {
    ifstream input;
    input.open(file_name);
    if (!input) {
    	throw runtime_error("Error opening input file");
    }

    string line;
    shared_ptr<Keyword> key = nullptr;
    unordered_map<string, shared_ptr < Keyword >> ::iterator loc_finder;

    while (getline(input, line)) {
        line = trim(line);
        if (line == "END OF INPUT")
            break;
        loc_finder = input_keys.find(line);
        if (loc_finder != input_keys.end()) {
            // this is a keyword
            key = loc_finder->second;
            // check that we haven't set it yet
            if (key->isValueSet())
                throw runtime_error("Duplicate keyword: " + key->getName());
            // read in the value(s)
            while (!key->isValueSet()){
                key->parse(input);
            }
        }
        else
            throw runtime_error("Unknown keyword: " + line);
    }

    // check that we haven't missed any required keywords
    for (auto key_val : input_keys) {
        key = key_val.second;
        if (key->isRequired() && !key->isValueSet())
            throw runtime_error("The required keyword " + key_val.first +
                    " was not found in the input file");
    }

}

// print out all keywords and their values
void InputReader::output_keywords() {
    
    cout << endl;
    cout << "********************************************************************************" << endl;
    cout << "**********                                                            **********" << endl;
    cout << "                        CONFIGURATION OF THIS CALCULATION                       " << endl;
    cout << endl;
    cout << setw(40) << left << "Keyword" << "Value" << endl;
    cout << setw(40) << left << "-------" << "-----" << endl;


    // let's make the output as easy to parse as possible
    // by writing out the keywords in alphabetical order
    vector<string> keys;
    for (auto key_val : input_keys)
        keys.push_back(key_val.first);
    sort(keys.begin(), keys.end());

    for (string key : keys) {
        // keyword name
        cout << setw(40) << left << key;

        // keyword value
        auto val = input_keys[key];
        // check if int keyword
        shared_ptr<IntKeyword> ik = dynamic_pointer_cast<IntKeyword> (val);
        if (ik) {
            cout << ik->getIntKeyword();
        }

        // check if string keyword
        shared_ptr<StringKeyword> sk = dynamic_pointer_cast<StringKeyword> (val);
        if (sk) {
            cout << sk->getStringKeyword();
        }

        // check if bool keyword
        shared_ptr<BoolKeyword> bk = dynamic_pointer_cast<BoolKeyword> (val);
        if (bk) {
            cout << bk->getBoolKeyword();
        }

        // check if double keyword
        shared_ptr<DoubleKeyword> dk = dynamic_pointer_cast<DoubleKeyword> (val);
        if (dk) {
            cout << dk->getDoubleKeyword();
        }

        // check if IntVecs keyword
        shared_ptr<IntVecsKeyword> iv =
                dynamic_pointer_cast<IntVecsKeyword> (val);
        if (iv) {
            vector<vector<int>> vecs = iv->getIntVecsKeyword();
            if (!(iv->isValueSet())) {
                cout << "Not defined" << endl;
                continue;
            }
            bool init = false; // special formatting for first line
            for (vector<int> vec : vecs) {
                if (!init)
                    init = true;
                else
                    cout << "                                        ";
                for (int i : vec) {
                    cout << i << "  ";
                }
                cout << endl;
            }
            continue;
        }

        // check if DoubleVecs keyword
        shared_ptr<DoubleVecsKeyword> dv =
                dynamic_pointer_cast<DoubleVecsKeyword> (val);
        if (dv) {
            vector<vector<double>> vecs = dv->getDoubleVecsKeyword();
            if (!(dv->isValueSet())) {
                cout << "Not defined";
                if (dv->hasDefault()) {
                    cout << ", use default:" << endl;
                    cout << "                                        ";
                } else {
                    cout << endl;
                    continue;
                }
            }
            bool init = false; // special formatting for first line
            for (vector<double> vec : vecs) {
                if (!init)
                    init = true;
                else
                    cout << "                                        ";
                for (double x : vec) {
                    cout << x << "  ";
                }
                cout << endl;
            }
            continue;
        }
        cout << endl;
    }
    
    cout << "**********                                                            **********" << endl;
    cout << "********************************************************************************" << endl;
}

// return value for string keyword given by name
string InputReader::getStringKeyword(string name){
    shared_ptr<Keyword> k = input_keys.at(name);
    shared_ptr<StringKeyword> sk =
            dynamic_pointer_cast<StringKeyword>(k);
    return sk->getStringKeyword();
}

// return value for int keyword given by name
int InputReader::getIntKeyword(string name){
    shared_ptr<Keyword> k = input_keys.at(name);
    shared_ptr<IntKeyword> ik = dynamic_pointer_cast<IntKeyword>(k);
    return ik->getIntKeyword();
}

// return value for bool keyword given by name
bool InputReader::getBoolKeyword(string name){
    shared_ptr<Keyword> k = input_keys.at(name);
    shared_ptr<BoolKeyword> bk = dynamic_pointer_cast<BoolKeyword>(k);
    return bk->getBoolKeyword();
}

// return value for double keyword given by name
double InputReader::getDoubleKeyword(string name){
    shared_ptr<Keyword> k = input_keys.at(name);
    shared_ptr<DoubleKeyword> dk =
            dynamic_pointer_cast<DoubleKeyword>(k);
    return dk->getDoubleKeyword();
}

// return value for IntVecs keyword given by name
vector<vector<int>> InputReader::getIntVecsKeyword(string name){
  shared_ptr<Keyword> k = input_keys.at(name);
  shared_ptr<IntVecsKeyword> iv = dynamic_pointer_cast<IntVecsKeyword>(k);
  return iv->getIntVecsKeyword();
}

// return value for DoubleVecs keyword given by name
vector<vector<double> > InputReader::getDoubleVecsKeyword(string name){
    shared_ptr<Keyword> k = input_keys.at(name);
    shared_ptr<DoubleVecsKeyword> dv = dynamic_pointer_cast<DoubleVecsKeyword>(k);
  return dv->getDoubleVecsKeyword();
}

// return a vector containing all keyword names
vector<string> InputReader::getNames() {
    vector<string> names;
    for (pair<string, shared_ptr<Keyword>> key : input_keys)
        names.push_back(key.first);
    return names;
}
