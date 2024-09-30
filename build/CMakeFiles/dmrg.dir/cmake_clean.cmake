file(REMOVE_RECURSE
  "libdmrg.a"
  "libdmrg.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/dmrg.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
