package require psfgen

topology ./top_all36_prot.rtf
pdbalias atom ILE CD1 CD
pdbalias residue HIS HSD

segment B {
    pdb my_7c8j.pdb
}
coordpdb my_7c8j.pdb B

guesscoord
writepsf "my_7c8j.psf"
writepdb "my_7c8j_processed.pdb"
exit
