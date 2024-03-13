import numpy as np

def Chain1BodyHamil( Nsites, t = 1., mu = 0., pbc = True ):
  """
  Returns the 1-body Hamiltonian for a 1d chain with Nsites sites,
  possibly periodic boundary conditions (pbc), hopping amplitude t
  and chemical potential mu.
  """
  h = -mu * np.identity(Nsites, dtype = np.float64)
  for i in range(1, Nsites):
    h[i-1, i  ] = -t
    h[i  , i-1] = -t
  if( pbc ):
    h[0 , -1] = -t
    h[-1,  0] = -t

  return h

def GetGSEn( Nsites, Nup, Ndo, t = 1., mu = 0., pbc = True ):
  """
  Get the ground state energy of a 1d chain with Nsites sites,
  nearest neighbor hopping t, chemial potential mu, optional
  periodic boundary conditions (pbc), and Nup electrons of spin
  up, Ndo electrons of spin down.
  """
  h            = Chain1BodyHamil( Nsites, t, mu, pbc )
  evals, evecs = np.linalg.eigh(h)

  E0 = np.sum( evals[:Nup] ) + np.sum( evals[:Ndo] )
  return E0

def PlotEnergConv_HalfFilling( Nsites_max = 150, t = 1., mu = 0. ):
  """
  Plot energy converegence for finite size chains, with and without
  periodic boundary conditions.
  """
  nsites    = np.arange( 2, Nsites_max, 2 )
  pbc_E0s   = np.zeros_like( nsites, dtype = np.float64 )
  nopbc_E0s = np.zeros_like( nsites, dtype = np.float64 )

  for (i, nsite) in enumerate(nsites):
    nup          = nsite // 2
    ndo          = nsite - nup
    pbc_E0s[i]   = GetGSEn( nsite, nup, ndo, t, mu, pbc = True  ) / float(nsite)
    nopbc_E0s[i] = GetGSEn( nsite, nup, ndo, t, mu, pbc = False ) / float(nsite)

  import matplotlib.pylab as plt

  fig = plt.figure()
  ax  = fig.add_subplot(111)

  ax.plot( nsites,   pbc_E0s, linestyle = '--', marker = 'o', label = 'PBC'    )
  ax.plot( nsites, nopbc_E0s, linestyle = '--', marker = 'x', label = 'No PBC' )
  ax.axhline( -4. / np.pi, color = 'grey', label = 'Thermo. Limit' )
  
  ax.set_xlabel(r"Chain Length")
  ax.set_ylabel(r"Ground state energy per site $[t]$")
  
  ax.tick_params( which = "both", direction = 'in' )

  ax.legend(loc = 'best')
  plt.grid( alpha = 0.5, color = 'grey' )
  
  fig.savefig("nonint_en_conv.pdf") 

if __name__ == "__main__":

  import sys
  nsites = int(sys.argv[1])

  nup = nsites // 2
  ndo = nsites - nup

  print( "Nsites: ", nsites, ", nup: ", nup, ", ndo: ", ndo)
  Epbc = GetGSEn( nsites, nup, ndo, pbc = True ) / float(nsites)
  Eobc = GetGSEn( nsites, nup, ndo, pbc = False ) / float(nsites)

  print( "GS energy per site ( w   pbcs ): ", Epbc )
  print( "GS energy per site ( w/o pbcs ): ", Eobc )

  file = open("Energy_per_site.dat", "a")  # write mode
  file.write('{}\t{}\t{}\n'.format(nsites,Eobc,Epbc))
  file.close()
  
  PlotEnergConv_HalfFilling( Nsites_max = 150 )
