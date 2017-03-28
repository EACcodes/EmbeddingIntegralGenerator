***,C
print,basis,orbitals

symmetry,nosym
geometry={
4
angstrom
C   3.594222       6.218747       8.509041
O   2.916401       6.422257       7.534980
H   3.528797       5.268722       9.100681
H   4.337350       6.964926       8.891556
}

basis=cc-pVTZ

hf
en_no_emb = energy

{matrop
LOAD,one_int,H0
READ,emb_int,FILE=emb_ints.molpro
ADD,H0,one_int,emb_int
SAVE,H0,1200.1,h0
}

hf
en_emb = energy

show en_no_emb
show en_emb
