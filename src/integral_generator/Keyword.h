/*
 * File:   Keyword.h
 * Authors: Francis Ricci and David Krisiloff
 *
 */

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>


#ifndef KEYWORD_H
#define	KEYWORD_H

using namespace std;

// an overarching class for the storage of keywords and their values, containing
// many sub-classes
class Keyword {
public:
    // parse keyword from input file stream
    virtual void parse(ifstream& input)=0;
    
    // return true if value has been set, false otherwise
    bool isValueSet() {return _finished_reading; } ;
    
    // return true if keyword input is required, false otherwise
    bool isRequired() {return _essential;} ;
   
    // return keyword name
    string getName() {return _name;};

protected:
    // All keywords contain a name and booleans for set value and required input
    // Must provide name and required status on construction
    Keyword(string name, bool required):
    _essential(required), _finished_reading(false), _name(name) {}
    bool _essential;
    bool _finished_reading = false;
    string _name;
};


/*
 * Integer keywords come in two types, those with good initial defaults
 * and those without (and are therefore required!)
 */
class IntKeyword : public Keyword{
public:
    // Int keywords consist of a keyword plus an int value
    // Must provide name, default value, and required status on construction
    IntKeyword(string name, int default_val, bool required);
    
    // Specific parsing for int keywords
    void parse(ifstream& input);
    
    // return the value of this int keyword
    int getIntKeyword();

protected:
    // value of the int keyword
    int _val;
};

// required int keywords have no default value
class RequiredIntKeyword : public IntKeyword{
public:
    // constructor requires only the keyword name
    RequiredIntKeyword(string name);
};

// default int keywords have a good default value
class DefaultIntKeyword : public IntKeyword{
public:
    // constructor requires the name of the keyword and its default value
    DefaultIntKeyword(string name, int default_value);
};

/*
 * String keywords 
 */
class StringKeyword : public Keyword{
public:
    // constructor requires the name of the keyword
    StringKeyword(string name, string default_val, bool required);
    
    // specific parsing for string keywords
    void parse(ifstream& input);
    
    // return the value of the string keyword
    string& getStringKeyword();
private:
    // value of the string keyword
    string _val ; 
};

/*
 * With the density fitting work we now have optional string keywords
 */
class DefaultStringKeyword : public StringKeyword{
public:
    DefaultStringKeyword(string name, string default_val);
};

class RequiredStringKeyword : public StringKeyword{
public:
    RequiredStringKeyword(string name);
};


/*
 * Boolean keywords are either required input or come with a default value
 * Some boolean keywords take no arguments and are provided as just flags
 */
class BoolKeyword : public Keyword{
public:
    // constructor requires the name of the keyword, its default value, and
    // whether or not it is a required input
    BoolKeyword(string name, bool default_val, bool required);
    
    // specific parsing for boolean keywords
    void parse(ifstream& input);
    
    // return value of the boolean keyword
    bool getBoolKeyword();
    
protected:
    // value of the boolean keyword
    bool _val;
};

// required bool keywords have no default value
class RequiredBoolKeyword : public BoolKeyword{
public:
    // constructor requires keyword name
    RequiredBoolKeyword(string name);
};

// default bool keywords have a good default value
class DefaultBoolKeyword : public BoolKeyword{
public:
    // constructor requires keyword name and default value
    DefaultBoolKeyword(string name, bool default_val);
};

// flag bool keywords take no arguments
class FlagBoolKeyword : public BoolKeyword{
public:
    // constructor requires keyword name
    FlagBoolKeyword(string name);
    
    // specific parsing for flag bool keywords
    void parse(ifstream& input);
};

// default double keywords have a good default value
class DoubleKeyword : public Keyword{
public: 
    // constructor requires keyword name and default value
    DoubleKeyword(string name, double default_val);
    
    // specific parsing for double keywords
    void parse(ifstream& input);
    
    // return the value of the double keyword
    double getDoubleKeyword();
    
protected:
    // value of the double keyword
    double _val;
};

// IntVecs keywords are specific to the storage of integer vectors
class IntVecsKeyword : public Keyword{
public:
    // constructor requires keyword name
    IntVecsKeyword(string name, bool required);

    // specific parsing for integer vectors
      void parse(ifstream& input);

    // returns a reference to the vector of integers
    vector<vector<int>>& getIntVecsKeyword();

protected:
    // value of the integer vector
    vector<vector<int>> _val;
};

// required vector of double vectors 
class RequiredIntVecsKeyword : public IntVecsKeyword{
public:
    // constructor requires also default
    RequiredIntVecsKeyword(string name);
};
 
// DoubleVecs keywords are specific to the storage of vectors
// containing vectors of doubles
class DoubleVecsKeyword : public Keyword{
public:
    // constructor requires keyword name
    DoubleVecsKeyword(string name, bool required);

    // specific parsing for double matrices
    void parse(ifstream& input);

    // returns a reference to the matrix of doubles
    vector<vector<double>>& getDoubleVecsKeyword();

    // return true if keyword input has a default, false otherwise
    bool hasDefault() {return _default;} ;

protected:
    // constructor requires keyword name and default value
    DoubleVecsKeyword(string name, vector<vector<double>> default_value);
    // value of matrix of doubles
    vector<vector<double>> _val;
    vector<vector<double>> _default_val;
    bool _default;
};

// default vector of double vectors 
class DefaultDoubleVecsKeyword : public DoubleVecsKeyword{
public:
    // constructor requires also default
    DefaultDoubleVecsKeyword(string name, vector<vector<double>> default_value);
};

// required vector of double vectors 
class RequiredDoubleVecsKeyword : public DoubleVecsKeyword{
public:
    // constructor requires also default
    RequiredDoubleVecsKeyword(string name);
};




#endif	/* KEYWORD_H */
