# Lattice DMRG 

This is a simple, yet complete, **parallel DMRG library** to solve interacting spin-S (Heisenberg) and fermions (Hubbard) models in 1D. The software exploits consevation of the Quantum Numbers (e.g. $S_z$ or $[N_\uparrow,N_\downarrow]$), high-performance objects/algorithms and distributed MPI framework to perform *infinite* and *finite* DMRG algorithms.  
 
The structure of this code is largely inspired by the simple-DMRG project: [GitHub](https://github.com/simple-dmrg/simple-dmrg) and [Zenodo](https://zenodo.org/record/1068359).

### Dependencies
The code is based on:  

- [X] [SciFortran](https://github.com/aamaricci/SciFortran)  

- [X] MPI 

  
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

#### <a name="milestone5"></a> Milestone 5 
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
- [X] Further development of *finite* DMRG algorithm.


#### Milestone 7 (final?)
- [X] Massive distributed MPI parallelization of the DMRG algorithm.   

So far different ideas have been be explored in the literature. Some focused on the parallelization of the DMRG algorithm itself, by distributing workload from different sites/MPS across nodes. Other have focused on accelarating the parallel execution of the tensor product for $H_{sb}$ for `sparse_H=T`, see for instance Ref.[[1]](#1). 

Here we take a slightly different approach. A simple profiling shows that more than 80% of the code executin time is spent in the solution of the eigenvalue problem for the superblock $H_{sb}|v\rangle = E|v\rangle$. 

As we have shown in [milestone 5](#milestone5), the super-block Hamiltonian $H_{sb}$ can be decomposed into different parts, each having the form of a suitable tensor product. This form is analog to the one assumed by Hubbard-like models with spin-conservation for which a massively parallel Exact Diagonalization algorithm is available, see Ref.[[2]](#2)-[[3]](#3).

Using this approach we constructed a simple yet efficient distributed MPI DMRG scheme as follows. We represent the generic superblock vector $|v\rangle \in {\cal H}_{sb}$ as a block matrix $V=diag\{V_q\}_{q=1,\dots,N_q}$ with $q$ the quantum number index and $N_q$ the total number of quantum number sectors contributing to the current superblock. Each block $V_q$ corresponds to a sub-matrix with $D_R(q)$ rows and $D_L(q)$ columns. 

In the current MPI impementation of the DMRG we distribute across the $M$ nodes a share of $Q=D_L(q)/M$ columns. As we will explain in the following, this enables to perform application of each term in $H_{sb}$ locally in the nodes $p$ memory with only minimal MPI communication overhead.  

For, we consider the function `vector_transpose_MPI` in `DMRG_GLOBAL`. This function exploits the MPI function `MPI_All2AllV` to perform the parallel transposition of a given matrix $A$ which is distributed column-wise. 

In addition we note that each term appearing in the expression for $H_{sb}$ corresponds to a given block with quantum number $q$, as by dictated by the symmetries of the problem: $H_{sb}=diag\{H_{sb}(q)\}$. Thus we can analyze the matrix-vector product (MVP) of each term with a given quantum number one-by-one:


* $1_L\otimes H_R$: when applied to the block $V_q$ this term involves multiplication of the (small) matrix $H_R$ rows with columns of $V_q$ with indices $j_q(p)$, with $p$ the node index. This multiplication is local in the memory of the node and does not required communication. 

* $H_L\otimes 1_R$: when applied to the block $V_q$ this term involves multiplication of the (small) matrix $H_L$ rows with the *rows* of $V_q$ with indices $i_q(p)$. As the rows are distributed among the nodes a global communication is required. To optimize this step we operate as follows: 
    * MPI-transpose the block $V_q\rightarrow V^T_q$, 
    * apply the matrix $H_L$ which now operates locally in the nodes memory
    * MPI-transpose back the result $(H_LV^T_q)^T\rightarrow H_LV$.

* $H_{LR} =  \sum_a \sum_{k} A_a(k)\otimes B_a(k)$: This term requires a more involved treatment. First we observe that each term in the double sum applies as:  

$$
[A_a(k)\otimes B_a(k)]|v_q\rangle = [A_a(k)\otimes B_a(k)]\cdot vec(V_q)
$$

where the operator $vec({\cdot})$ takes the vector on input and transforms it into a matrix columns-wise. By a well-known properties of the tensor products the last term is equal to:

$$
vec{(B_a(q,k)\cdot V_k \cdot A_a(k,q)^T)}
$$

which is an apparently involved double matrix-matrix product (MMP). However, as indicated in the last expression, because of the block structure of the operators $A_a$ and $B_a$ (indicated by the index $k$) some restrictions applies to this product: only the *off-diagonal* block components of $A_a$ and $B_a$ which ensures the final result of the MMP contributes to the specific $q$ quantum number are possible. 

Recalling that each $V_k$ is distributed column-wise we perform the first product in parallel as for $1_L\otimes H_R$: $C_a(q,k)= B_a(q,k)\cdot V_k$. The resulting dense  matrix $C_a(q,k)$ is  distributed column-wise by construction. 
We are then left with a final MMP:  

$$
vec{(C_a(q,k)\cdot A_a(k,q)^T )} = 
vec{( [A_a(q,k)\cdot C^T_a(k,q)]^T )}
$$

which, as for the term $H_L\otimes 1_R$ requires using MPI-transpose two times: i) to get $C^T_a(k,q)$ and ii) to transpose the final result to be accumulated in the outgoing vector. 

The overall cost of communication is minimized relying exclusively on the `MPI_All2All`-type of communication, thus unlocking massively parallel scaling of the MVP at the heart of the super-block diagonalization. 

  

### Results
Here are some results for the Heisenberg model:  
$H = J_z \sum_{i} S_z(i)S_z(i+1) + J_{xy} \sum_{i,a=x,y} S_a(i)S_a(i+1)$

![plot](https://github.com/aamaricci/Lattice_DMRG/blob/main/.plot/figs.png)

In the top panels we show the groundstate energy $E(j)$ and the entanglement entropy of the left block sites $S_L(j) = -{\mathrm Tr}[{\rho_L}\ln{\rho_L}]$ as a function the chain sites in a Spin 1/2 chain of 300 sites. The nearest-neighbor spin-spin correlation $\langle \vec{S}(j)\cdot \vec{S}(j+1)\rangle$ for the same system is reported in the bottom-left panel. Compare this result to the fig.6 of [[4]](#4).  
Finally, in the bottom-right panel we report the spatial distribution of the local magnetization for a Spin 1 chain with open boundary conditions, showing Spin 1/2 edge modes.


Results for the Hubbard model $H = -t\sum \sum_{i\sigma} c^\dagger_{i\sigma} c_{i+1\sigma} + 
U \sum_{i} n_{i\uparrow} n_{i\downarrow}$
are available here: 

![plot](https://github.com/QcmPlab/Lattice_DMRG/blob/main/.plot/figH.png)

In the top-left panel we compare the energy per site $E(j)$ with respect to the exact solution $E_0=-4t/\pi$ and the exact numberical solution for the non-interacting case at one electron per site $\langle N\rangle=1$. We used $M=40$ states to reach a satisfactory convergence of the ground state energy. In the top-right panel we report the occupation profile per spin. In the bottom-left panel we show the evolution of the entanglement entropy per site $S(j)$, while the last panel bottom-right displayes the profile of the double occupation.  


![gif](https://github.com/aamaricci/Lattice_DMRG/blob/main/.plot/DMRG_record.gif)




### Known issues
There are a number of known issues with this code which, mostly for time reasons, we did not solve completely. Please report to any of the authors.
    
    
    
### References
<a id="1">[1]</a> 
[Large-Scale Implementation of the Density Matrix Renormalization Group Algorithm](https://iris.sissa.it/handle/20.500.11767/68070), J.Vance, MHPC Thesis (2017). 

<a id="2">[2]</a> 
[EDIpack: A parallel exact diagonalization package for quantum impurity problems](https://doi.org/10.1016/j.cpc.2021.108261), A.Amaricci *et al.*, Computer Physics Communications, **273**, 108261 (April 2022).

<a id="3">[3]</a> 
[Massively parallel
exact diagonalization of strongly correlated systems](https://juser.fz-juelich.de/record/136300/files/Dolfen_Diplom.pdf), A.Dolfen, Master Thesis, IFF Julich (2006).

<a id="4">[4]</a> 
[Density-matrix algorithms for quantum renormalization groups](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.10345), S.White, Phys. Rev. B **48**, 10345 (1993).


### Info
For any information or suggestion contact the author OR create an issue in this repo.


--


The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


