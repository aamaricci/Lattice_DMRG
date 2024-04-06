import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

plt.rcParams["font.family"] = "Times"
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 28})
plt.rcParams['figure.figsize'] = [21, 14]




# ##################################################################
# ##################################################################
# HEISENBERG:
# ##################################################################
# ##################################################################
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, width_ratios=[1,1])


msize=12
ls = 'solid'
lw = 2

#Panel 1: Spin 1/2 E vs j L=300 M=40 #Exact Energy= -0.448949
##################################################################
ax1.text(150, -0.41, r'Energy per site Spin-$\frac{1}{2}$',fontsize=28,ha="center")
ax1.axis([0,300,-0.45,-0.4])               
ax1.set_xticks(np.arange(0,300, step=50))
ax1.set_xlabel(r'$j$')
ax1.set_ylabel(r'$E(j)$')
data = pd.read_csv('HeisenbergSU2.3/energyVSsys.length_L300_M40_iDMRG.dmrg', sep='\s+', header=None);x1=data[0];y1=data[1]
ax1.hlines(-0.448949,xmin=0,xmax=300,ls=':',lw=lw,colors='grey',alpha=0.8)
ax1.plot(2*x1,  y1,  lw=lw, ls=ls,  marker='o', mfc='red', ms=msize, c='blue')


#Panel 2: Spin 1/2  <S vs j L=300 M=40
##################################################################
ax2.text(75,0.4, r'Entropy per site Spin-$\frac{1}{2}$',fontsize=28,ha="center")
ax2.axis([0,150,0.3,1.15])               
ax2.set_xticks(np.arange(0,150, step=20))
ax2.set_xlabel(r'$j$')
ax2.set_ylabel(r'$S(j)$')
data = pd.read_csv('HeisenbergSU2.3/SentropyVSsys.length_L300_M40_iDMRG.dmrg', sep='\s+', header=None);x1=data[0];y1=data[1]
ax2.plot(x1,  y1,  lw=lw, ls=ls, marker='o', mfc='red', ms=msize, c="blue")


#Panel 3: Spin 1/2 <S_iS_{i+1} vs j L=100 M=40
##################################################################
ax3.text(75,-0.3225, r'Spin-Spin correlation Spin-$\frac{1}{2}$',fontsize=28,ha="center")
ax3.axis([0,150,-0.7,-0.25])               
ax3.set_xticks(np.arange(0,150, step=20))
ax3.set_xlabel(r'$j$')
ax3.set_ylabel(r'$\langle \vec{S}(j)\cdot \vec{S}({j+1}) \rangle$')
data = pd.read_csv('HeisenbergSU2.3/SiSjVSjL300_M40_User.dmrg', sep='\s+', header=None);x1=data[0];y1=data[1]
ax3.plot(x1,  y1,  lw=lw, ls=ls,  marker='o', mfc='red', ms=msize, c='blue')



#Panel 4: Spin 1 <S_z(i) vs j L=100 M=40
##################################################################
ax4.text(25,0.45, r'Local magnetization Spin-1',fontsize=28,ha="center")
ax4.axis([0,50,-0.4,0.6])               
ax4.set_xticks(np.arange(0,50, step=10))
ax4.set_xlabel(r'$j$')
ax4.set_ylabel(r'$\langle S_z(j)\rangle$')
data = pd.read_csv('HeisenbergSU2.3/SzVSj_L100_Err1.000000E-05_Spin1.dmrg', sep='\s+', header=None);x1=data[0];y1=data[1]
ax4.plot(x1,  y1,  lw=lw, ls=ls,  marker='o', mfc='red', ms=msize, c='blue')


plt.savefig("figs.pdf", bbox_inches='tight')
plt.savefig("figs.png", bbox_inches='tight')
plt.clf()
plt.cla()
plt.close()





# ##################################################################
# ##################################################################
# HUBBARD
# ##################################################################
# ##################################################################

fig = plt.figure()
ax  = plt.gca()
ax.text(35, -1.125, r'Energy per site Hubbard $U=0$, $\langle N\rangle=1$',fontsize=28,ha="center")
ax.axis([0,100,-1.28,-1.11])               
ax.set_xticks(np.arange(0,100, step=10))
ax.set_xlabel(r'$j$')
ax.set_ylabel(r'$E(j)$')
data = pd.read_csv('Hubbard/TBexact/eVSj.obc', sep='\s+', header=None);x0=data[0];y0=data[1]
data = pd.read_csv('Hubbard/sparse_True/energyVSleft.length_L100_M100_iDMRG.dmrg', sep='\s+', header=None);x1=data[0];y1=data[1]
data = pd.read_csv('Hubbard/sparse_False/energyVSleft.length_L100_M100_iDMRG.dmrg', sep='\s+', header=None);x2=data[0];y2=data[1]

msize=12
ls = 'solid'
lw = 2

ax.hlines(-4/np.pi,xmin=0,xmax=100,ls=':',lw=lw,colors='grey',alpha=1.0,label=r"$-4t/\pi$")
ax.plot(2*x0,  y0,  lw=lw, ls=ls, c='black', label=r"Tight Binding")
ax.plot(2*x2,  y2, lw=lw, ls="dashed",  marker='x', mfc='green', ms=14, c='blue', label=r"Direct H*v")
ax.plot(2*x1,  y1,  lw=lw, ls='None', marker='o', mfc='red', ms=10,  alpha=0.5,label=r"Sparse H*v")
ax.legend(loc='upper right', fontsize = '24', handlelength=1.0, framealpha=0,ncol=2);
plt.savefig("figH.pdf", bbox_inches='tight')
plt.savefig("figH.png", bbox_inches='tight')
plt.clf()
plt.cla()
plt.close()
