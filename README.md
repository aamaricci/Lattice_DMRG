# Lattice DMRG

*beta version*  
*This code mostly serves as development platform for a Quantum Impurity Solver algorithm.*


A simple DMRG code to solve Heisenberg (spin-S) and Hubbard models in 1d.  
The software exploits consevation of the Quantum Numbers tuple ($S_z$ or $[N_\uparrow,N_\downarrow]$) to perform *infinite* and *finite* DMRG algorithms.  

 
The structure of this code is largely inspired by the simple-DMRG project: [GitHub](https://github.com/simple-dmrg/simple-dmrg) and [Zenodo](https://zenodo.org/record/1068359).

### Dependencies

The code is based on:  

- [X] [SciFortran](https://github.com/aamaricci/SciFortran)  

- [ ] MPI 

  
### Installation

Clone the repo  
Open and setup the Makefile with your favorite editor (hint: `emacs`)  
Compile  

```
git clone https://github.com/aamaricci/Lattice_DMRG 
```
```
cd Lattice_DMRG
emacs Makefile
```
```
make
```


### Info
For any information or suggestion contact the authors:   
adriano DOT amaricci AT gmail DOT com  
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
- [x] Wrap the SuperBlock construction into a dedicated module. Encapsulate and generalize the code structure. 
- [x] Merge with Spin DMRG code. 


#### Milestone 4
- [x] Add layer to save rotation/truncation matrices
- [x] Develop *finite* DMRG algorithm
- [x] Measure local observables 
- [ ] Measure entanglement
- [ ] Implement a better strategy for the matrix-vector product H_sb*V, using tensor product structure of H_sb. 

#### Milestone 5
- [ ] Implement parallel tensor product. 

#### Future developemnts
- [ ] Development of a iDMRG algorithm for fermions exploiting the spin separability of the Hamiltonian. See [https://doi.org/10.1016/j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261)


--

***LICENSE***  
Copyright 2023- (C) Adriano Amaricci, Carlos Mejuto Zaera, Ricardo Lopes

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


