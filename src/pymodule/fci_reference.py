# Minimal FCI for H2 in minimal basis (sto-3g)
# Reference implementation for validation of VB-CI
import numpy as np
from itertools import combinations

def build_fci_h2_hamiltonian(h_core, eri, n_electrons=2):
    """
    Build the FCI Hamiltonian for H2 in minimal basis (2 orbitals, 2 electrons, spin-restricted).
    h_core: (2,2) core Hamiltonian in MO basis
    eri: (2,2,2,2) two-electron integrals in MO basis (chemist's notation)
    Returns: (Ndet, Ndet) Hamiltonian matrix, list of determinants (as tuples of occupied orbital indices)
    """
    # All possible doubly-occupied determinants (alpha/beta pairs)
    # For 2 electrons in 2 spatial orbitals, there are 3 singlet determinants:
    # |1a 1b>, |2a 2b>, (1/sqrt(2))(|1a 2b> + |2a 1b>)
    dets = [(0,0), (1,1), (0,1)]
    N = len(dets)
    H = np.zeros((N,N))
    # Diagonal elements
    for i, d in enumerate(dets):
        occ = list(set(d))
        H[i,i] = 2 * h_core[occ[0], occ[0]] if len(occ)==1 else sum(h_core[o,o] for o in occ)
        H[i,i] += eri[occ[0],occ[0],occ[0],occ[0]] if len(occ)==1 else sum(eri[o,o,p,p] for o,p in combinations(occ,2))
    # Off-diagonal (only (0,0)<->(0,1) and (1,1)<->(0,1) coupled by single excitation)
    H[0,2] = H[2,0] = np.sqrt(2) * h_core[0,1]
    H[1,2] = H[2,1] = np.sqrt(2) * h_core[1,0]
    return H, dets

def solve_fci_h2(h_core, eri):
    H, dets = build_fci_h2_hamiltonian(h_core, eri)
    E, C = np.linalg.eigh(H)
    idx = np.argmin(E)
    return E[idx], C[:,idx], H, dets
