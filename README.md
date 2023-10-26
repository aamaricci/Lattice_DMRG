# Lattice DMRG
*beta version*  

This is simple, yet complete, DMRG code to solve Heisenberg (spin-S) and Hubbard models, prevalently in 1d. The software exploits consevation of the Quantum Numbers ($S_z$ or $[N_\uparrow,N_\downarrow]$) to perform *infinite* and *finite* DMRG algorithms.  

 
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
- [x] Measure entanglement
- [x] Measure nearest neighbor correlations

#### Milestone 5
- [ ] Implement a better strategy for the SuperBlock matrix-vector product $H_{sb}|\psi\rangle$, exploiting the tensor product structure of $H_{sb}= H_L\otimes 1_R + 1_L\otimes H_R + H_{LR}$. 
- [ ] Implement parallel tensor product, see [Large-Scale Implementation of the Density Matrix Renormalization Group Algorithm](https://iris.sissa.it/handle/20.500.11767/68070). 

#### Future developemnts
- [ ] Development of a iDMRG algorithm for fermions exploiting the spin separability of the Hamiltonian. See [https://doi.org/10.1016/j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261)


### Results
Here are some results for the Heisenberg model:  
$H = J_z \sum_{i} S_z(i)S_z(i+1) + J_{xy} \sum_{i,a=x,y} S_a(i)S_a(i+1)$

![plot](https://github.com/QcmPlab/Lattice_DMRG/blob/main/plot/figs.png)

In the top panels we show the groundstate energy $E(j)$ and the entanglement entropy of the left block sites $S_L(j) = -{\mathrm Tr}[{\rho_L}\ln{\rho_L}]$ as a function the chain sites in a Spin 1/2 chain of 300 sites. The nearest-neighbor spin-spin correlation $\langle \vec{S}(j)\cdot \vec{S}(j+1)\rangle$ for the same system is reported in the bottom-left panel. Compare this result to the fig.6 of [Density-matrix algorithms for quantum renormalization groups](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.10345) (S.White, Phys. Rev. B 48, 10345).  
Finally, in the bottom-right panel we report the spatial distribution of the local magnetization for a Spin 1 chain with open boundary conditions, showing Spin 1/2 edge modes.

![gif](https://github.com/QcmPlab/Lattice_DMRG/blob/main/plot/DMRG_record.gif)

### Known issues
There are a number of known issues with this code which, mostly for time reasons, we did not solve completely. Please report to any of the authors.
    
    
--


***LICENSE***  
Copyright 2023- (C) Adriano Amaricci, Carlos Mejuto Zaera, Ricardo Lopes and Massimo Capone

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


