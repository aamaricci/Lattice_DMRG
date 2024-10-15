MODULE DMRG
  !GENERAL
  USE INPUT_VARS
  USE AUX_FUNCS
  !AUX MODULES
  USE MATRIX_SPARSE
  USE MATRIX_BLOCKS
  USE TUPLE_BASIS
  USE LIST_OPERATORS
  USE LIST_SECTORS
  !DMRG:
  USE SITES
  USE BLOCKS
<<<<<<< HEAD
<<<<<<< HEAD:src/DMRG.f90
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
  USE DMRG_SUPERBLOCK
  USE DMRG_MAIN
  USE DMRG_MEASURE
  !
  !USE VARS_GLOBAL!< uncomment this to compile kron_ tests
<<<<<<< HEAD
=======
  USE SYSTEM
  USE MEASURE
>>>>>>> f63915b (Testing the code.):DMRG.f90
=======
>>>>>>> 7e90d6a (Updating Cmake library construction)
  implicit none  
END MODULE DMRG
