***,C
print,basis,orbitals

symmetry,nosym
geometry={
1
angstrom
C 4.00 5.00 5.00
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
