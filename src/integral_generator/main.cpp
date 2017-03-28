/*
 * File: main.cpp
 * Authors: Francis Ricci, David Krisiloff, Johannes M Dieterich
 *
 */

#include <iostream>
#include <string>
#include <stdexcept>

#include "integral_generator/InputReader.h"
#include "integral_generator/Keyword.h"
#include "integral_generator/ConvertKeys.h"
#include "integral_generator/Timer.h"

#include "integral_generator/BasisData.h"
#include "integral_generator/EmbeddingData.h"
#include "integral_generator/TigerStructs.h"

#ifdef _OPENMP
#include <omp.h>
#endif 

using namespace std;

int main(int argc, char** argv) {
    
        if(argc == 1){
        cout << "ERROR: EMBEDDING_INTEGRAL_GENERATOR must be called with an input file as first argument." << endl;
        return 1;
    }
        
    Timer total_time;
        
        cout << endl;
        cout << "********************************************************************************" << endl;
        cout << "********************************************************************************" << endl;
        cout << "**********                                                            **********" << endl;
        cout << "                        EMBEDDING INTEGRAL GENERATOR                            " << endl;
        cout << "**********                                                            **********" << endl;
        cout << "********************************************************************************" << endl;
        cout << "********************************************************************************" << endl;

    ConvertKeys global_vars;
    try {
        // read input file
        InputReader setup;
        setup.read_tiguar_input(argv[1]);
        setup.output_keywords();

        // convert keywords into global variables
        global_vars.keys_to_globals(setup);
    }
    catch (runtime_error &e){
        cerr << "FATAL ERROR! " << e.what() << endl;
        return tiger::INPUT_FAILURE;
    }

#ifdef _OPENMP
    // tell OMP how many threads to use
    omp_set_num_threads(global_vars.get_global_int("numThreads"));
    cout << endl << "Running with OpenMP using " << global_vars.get_global_int("numThreads");
    cout << " threads" << endl;
#endif

    BasisData *basis;
    try{
        basis = new EmbeddingBasis(global_vars);
    }
    catch (runtime_error &e){
        delete basis;
        cerr << "FATAL ERROR! " << e.what() << endl;
        return tiger::BASIS_FAILURE;
    }


//    BasisData *basis;
//    try{
//        basis = new ErkaleBasis(global_vars, tiger::full);
//    }
//    catch (runtime_error &e){
//        delete basis;
//        cerr << "FATAL ERROR! " << e.what() << endl;
//        return tiger::BASIS_FAILURE;
//    }

    EmbeddingData *embpot;
    try{
        embpot = new EmbeddingData(global_vars, basis);
    }
    catch (runtime_error &e){
        delete basis;
        delete embpot;
        cerr << "FATAL ERROR! " << e.what() << endl;
        return tiger::EMB_INT_FAILURE;
    }

    try{
        embpot->compute_integrals();
    }
    catch (runtime_error &e){
        delete basis;
        delete embpot;
        cerr << "FATAL ERROR! " << e.what() << endl;
        return tiger::EMB_INT_FAILURE;
    }

    delete embpot;
    delete basis;
    // Just before exiting print total execution time
    total_time.print_elapsed_time("TOTAL EXECUTION TIME ");
    return tiger::SUCCESS;
}

