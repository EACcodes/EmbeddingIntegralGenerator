# EmbeddingIntegralGenerator

Generator of embedding integrals in different formats to allow for reading
with different quatum chemistry codes like TigerCI, Molpro, Molcas, etc.

Currently available (To be continued):
* Molpro
* TigerCI
* GAMESS
* MOLCAS
* QChem

Needs basis set specification, molecular structure, embedding potential, 
and grid information. 

## Compilation

Cmake based build system. Intel compilers can be replaced with gnu compilers but not tested.

1. Configure with

        CC=icc CXX=icpc cmake <path to src> -DCMAKE_BUILD_TYPE=<Debug or Release> -DOMP=<On or Off>

    For production calculations use -DCMAKE_BUILD_TYPE=Release and -DOMP=On

2. Compile with

        make -j <number of threads>

## Usage

To run the code call it with an input file as parameter

    emb_ints.exe <input file> > <output file>

Currently, two output formats are possible (TIGERCI and MOLPRO, see below for details).
For a list of keywords and examples see below.

## Examples

* Example input files are included in the folder "examples". 
* There are two examples each for every output format.
* In the case of Molpro, examples for interfacing with a HF calculation are included 
  and the files "reference\_energies\_molcas" contain HF energies obtained with and without 
  the potential from Molcas for comparison.
* The embedding potential applied in these calculations was generated from random numbers.

## Tests

To be added. So far only files with reference energies from Molcas have been included with the
example for the Molpro interface.

## List of keywords

* GRID DIMENSIONS
    * Three integers on the following line to specify grid dimension in
      x, y, z direction
    * Required

* LATTICE VECTORS 
    * The lattice vectors in x, y, z, direction on the next three lines
    * Required

* SHIFT VECTOR
    * Specifies the shift vector on the next line
    * Default: 0.0  0.0  0.0
    * Use with care, not well tested

* XYZ FILE
    * xyz file containing geometry information
    * Required

* BASIS FILE 
    * file with basis set in Molcas format 
    * all-electron nuclear charge has to be specified even when using ECPs (old format)
    * include the basis set only (remove the ECP block)
    * Required

* EMBEDDING POTENTIAL FILE
    * file containing the embedding potential on the specified grid
    * Required

* OUTPUT FORMAT
    * string specifying the desired output format for the embedding integrals
    * available: MOLPRO, TIGERCI, GAMESS, MOLCAS

* AO INTEGRAL THRESHOLD (not working yet)
    * Planned to be threshold for printing integrals
    * Default: 1.0e-12

* CARTESIAN
    * Flag switching to cartesian instead of spherical basis functions
    * Not fully supported use with care. Will be improved in  future.

* EMBEDDING INTS REAL
    * Flag switching to old (depreciated) version of the integration scheme
      purely in real space
    
* NUM THREADS
    * Number of threads to use
    * Default: 1

