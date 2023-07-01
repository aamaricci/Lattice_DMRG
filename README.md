# Lattice DMRG

*beta version*

This is a simple DMRG code to solve interacting electrons problem on a lattice (1d so far). It mostly serve as a milestone in the development of a Quantum Impurity solver for Dynamical Mean-Field Theory. 

The algorithm includes consevation of Quantum Number $N_\up$, $N_\dw$.  
 

The structure of this code is largely based on the simple Spin_DMRG and extend it to the case of interacting electrons. It is possible in the future the two codes will merge.  

### Dependencies

The code is based on:  

* SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  

* MPI [ in a future release ]

  
### Installation

Installation is not fully supported so far. In a future release a CMake installation will be made available.    

For the time being Clone the repo and DIY


For any information or suggestion contact the author as:  
adriano DOT amaricci @ gmail DOT com

or 

create an issue in this repo.

### DEVELOPMENT

#### Milestone 1
- [x] Extend Spin iDMRG code to Electronic case. Solve the fermion sign problem.
- [x] Develop test code with recursive block construction and digonalization.
- [x] Develop code without Quantum Number (QN) conservationn. 
- [x] Test code against ED solution and exact behavior for the non-interacting regime.

#### Milestone 2
- [x] Develop the infrastructure to deal with QN conservation. 
- [x] Develop the construction of the states of a given QN sector. 
- [x] Code recursive blocking withing a specific QN sector. Compare with ED. 
- [x] Develop full iDMRG algorithm.
- [x] Test code against iDMRG without QN and ED.

#### Milestone 3
- [ ] Wrap the SuperBlock construction into a dedicated module. 
- [ ] Implement a better strategy for the matrix-vector product H_sb*V, using tensor product structure of  H_sb. 
- [ ] Implement parallel tensor product. 

#### Milestone 4
- [ ] Development of a iDMRG algorithm exploiting the $\up$,$\dw$ separability of the Hamiltonian. See [https://doi.org/10.1016/j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261)
- [ ] Include QN conservation following the results of **Milestone 2**
- [ ] Develop massively parallel scaling.


--

***LICENSE***  
Copyright 2023- (C) Adriano Amaricci, Carlos Mejuto Zaera

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


