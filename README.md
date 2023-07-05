# Hubbard DMRG

*beta version*

A simple DMRG code to solve the Hubbard model in 1d.  
*This code mostly serves as a milestone in the development of a Quantum Impurity solver for Dynamical Mean-Field Theory.*

The algorithm includes consevation of Quantum Number $N_\uparrow$, $N_\downarrow$ and perform *infinite* DMRG algorithm, although   
 

The structure of this code is based on the simple Spin_DMRG and extends that to the case of interacting electrons.   
Both are largely inspired by simple DMRG code: [GitHub](https://github.com/simple-dmrg/simple-dmrg) and [Zenodo](https://zenodo.org/record/1068359).

### Dependencies

The code is based on:  

- [X] SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  

* [ ] MPI 

  
### Installation

Clone the repo  
Setup the Makefile  
Compile


### Info
For any information or suggestion contact the authors: adriano DOT amaricci AT gmail DOT com
and/or
cmejutoz AT sissa DOT it 

OR 

create an issue in this repo.


### DEVELOPMENT

#### Milestone 1
- [x] Extend Spin iDMRG code to Electronic case from the Spin problem. Solve the fermion sign problem.
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
- [ ] Add layer to save rotation/truncation matrices
- [ ] Evaluates observables and entanglement
- [ ] Wrap the SuperBlock construction into a dedicated module. 
- [ ] Develop *finite* DMRG algorithm
- [ ] Implement a better strategy for the matrix-vector product H_sb*V, using tensor product structure of  H_sb. 
- [ ] Implement parallel tensor product. 

#### Future developemnts
- [ ] Development of a iDMRG algorithm exploiting the $\up$,$\dw$ separability of the Hamiltonian. See [https://doi.org/10.1016/j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261)


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


