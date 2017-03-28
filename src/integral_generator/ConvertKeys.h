/* 
 * File:   ConvertKeys.h
 * Author: Francis Ricci
 *
 * Created on June 10, 2014, 11:00 AM
 */

#ifndef CONVERTKEYS_H
#define	CONVERTKEYS_H

#include <string>

#include "InputReader.h"

using namespace std;

class ConvertKeys {
public:
    // default constructor
    ConvertKeys(){};
    
    // convert the keys in key_map into variable name-value pairs,
    // store in an unordered map
    void keys_to_globals(InputReader& key_map);
    string get_my_string(string name) const;
    int get_global_int(string name) const;
    bool get_global_bool(string name) const;
    string get_global_string(string name) const;
    double get_global_double(string name) const;
    vector<vector<int>> get_global_int_vecs(string name) const;
    vector<vector<double>> get_global_double_vecs(string name) const;

    // function to add global variables that old TigerCI needs
    void add_global_int(string name, int val);
    void add_global_bool(string name, bool val);
    
private:
    // check input parameters
    void check_input();
    unordered_map<string, string> my_strings;
};
#endif	/* CONVERTKEYS_H */
