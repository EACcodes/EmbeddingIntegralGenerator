/*
 * File:    InputReader.h
 * Authors: Francis Ricci and David Krisiloff
 *
 */

#ifndef INPUTREADER_H
#define	INPUTREADER_H

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>

#include "Keyword.h"

using namespace std;

// class for reading inputs from a file into an unordered map of keys
class InputReader {
public:
    //constructor -> initializes unordered map
    InputReader();
   
    // Read and store input parameters from string file_name
    void read_tiguar_input(string file_name);

    // print out all keywords and their values
    void output_keywords();
    
    // return value for string keyword given by name
    shared_ptr<Keyword> getKeyword(string name){return input_keys.at(name);};
    
    // return value for string keyword given by name
    string getStringKeyword(string name); 
    
    // return value for int keyword given by name
    int    getIntKeyword(string name);
    
    // return value for bool keyword given by name
    bool   getBoolKeyword(string name);
    
    // return value for double keyword given by name
    double getDoubleKeyword(string name);
    
    // return value for IntVecs keyword given by name
    vector<vector<int>> getIntVecsKeyword(string name);
  
    // return value for DoubleVecs keyword given by name
    vector<vector<double> > getDoubleVecsKeyword(string name);
    
    // return a vector containing all keyword names
    vector<string> getNames();
    
private:
    // unordered map for the storage of the input keywords
    unordered_map<string, shared_ptr<Keyword>> input_keys;
};

#endif	/* INPUTREADER_H */

