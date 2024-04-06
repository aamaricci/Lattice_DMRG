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
- [x] Implement a better strategy for the SuperBlock matrix-vector product $H_{sb}|\psi\rangle$, exploiting the tensor product structure of $H_{sb}= H_L\otimes 1_R + 1_L\otimes H_R + H_{LR}$. 

With this update we introduce two distinct method to evaluate and diagonalize the SuperBlock (SB) Hamiltonian $H_{sb}$, according to 
the value of the parameter `sparse_H=T,F`.  
If `T` the matrix $H_{sb}$ is evaluated as a sparse matrix performing the kroenecker products above.  
If `F` each non-trivial part of the the matrix $H_{sb}$, i.e. $H_L$, $H_R$ and $H_{LR}$ are directly applied to the SB vector $|\psi\rangle$.To improve execution each term is decomposed into its  blocks structure corresponding to conserved quantum numbers $\vec{Q}=[Q_1,\dots,Q_N]$. In addition, exploiting the fact that each SB vector can be decomposed into a matrix form of Left and Right states, $|\psi\rangle = \sum_{lr} c_{lr} |l\rangle|r\rangle$, we can operate directly as follows:  
$[H_L\otimes 1_R] \circ |\psi\rangle \rightarrow H_L \circ \sum_l c_{lr}|l\rangle$  
$[1_L\otimes H_R] \circ |\psi\rangle \rightarrow H_R \circ \sum_r c_{lr}|r\rangle$   
and  
$H_{LR} \circ |\psi\rangle \rightarrow  \sum_p \sum_{Q_i} A_p(Q_i)\otimes B_p(Q_i) |\psi\rangle$  
where the matrices $A_p$, $B_p$ are the connecting operators (i.e. $c_{i\sigma}$ or $S_a(i)$) restricted between two blocks with quantum numbers $Q_i$ and $Q_j=Q_i-1_{\sigma/a}$.   


#### Milestone 6
- [ ] Implement parallel strategies. Different ideas can be followed:  
i) parallel execution of the tensor product for $H_{sb}$ for `sparse_H=T`, see [Large-Scale Implementation of the Density Matrix Renormalization Group Algorithm](https://iris.sissa.it/handle/20.500.11767/68070).  
ii) parallel execution of the little matrix-vector and matrix-matrix products for `sparse_H=F`. We can also profit of the matrix structure of the SB vector as done in massively parallel Exact Diagonalization algorithm, see  [EDIpack: A parallel exact diagonalization package for quantum impurity problems](https://doi.org/10.1016/j.cpc.2021.108261).   
iii) distributed storage of the vectors and sparse matrices. 




### Results
Here are some results for the Heisenberg model:  
$H = J_z \sum_{i} S_z(i)S_z(i+1) + J_{xy} \sum_{i,a=x,y} S_a(i)S_a(i+1)$

![plot](https://github.com/QcmPlab/Lattice_DMRG/blob/main/.plot/figs.png)

In the top panels we show the groundstate energy $E(j)$ and the entanglement entropy of the left block sites $S_L(j) = -{\mathrm Tr}[{\rho_L}\ln{\rho_L}]$ as a function the chain sites in a Spin 1/2 chain of 300 sites. The nearest-neighbor spin-spin correlation $\langle \vec{S}(j)\cdot \vec{S}(j+1)\rangle$ for the same system is reported in the bottom-left panel. Compare this result to the fig.6 of [Density-matrix algorithms for quantum renormalization groups](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.10345) (S.White, Phys. Rev. B 48, 10345).  
Finally, in the bottom-right panel we report the spatial distribution of the local magnetization for a Spin 1 chain with open boundary conditions, showing Spin 1/2 edge modes.

![gif](https://github.com/QcmPlab/Lattice_DMRG/blob/main/.plot/DMRG_record.gif)



Results for the Hubbard model $H = -t\sum \sum_{i\sigma} c^\dagger_{i\sigma} c_{i+1\sigma} + 
U \sum_{i} n_{i\uparrow} n_{i\downarrow}$
are available here: 

![plot](https://github.com/QcmPlab/Lattice_DMRG/blob/main/.plot/figH.png)

For sake of simplicity we compare the energy per site of the non-interacting case $U=0$ at one electron per site $\langle N\rangle=1$ obtained from the two available Matrix-Vector methods (see above Milestone 5) against the exact Tight Binding and the exact solution $E_0 = -4t/\pi$. We used $M=100$ states which, as evident from the plot, are not enough to observe a decent convergence to the exact result but demonstrate the correct behavior of the DMRG algorithm. 


### Known issues
There are a number of known issues with this code which, mostly for time reasons, we did not solve completely. Please report to any of the authors.
    
    
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


