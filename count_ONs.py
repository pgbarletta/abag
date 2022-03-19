import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
HBondAtom = namedtuple('HBondAtom', ['chainID', 'chain_type',
                       'CDR', 'resSeq', 'resname', 'index', 'serial', 'element', 'is_sidechain'])
HBond = namedtuple('HBond', ['donor', 'acceptor'])
PolarCount = namedtuple('PolarCount', ['cdr_SC', 'cdr_BB', 'epi_SC', 'epi_BB'])
source_location = Path().resolve()
sys.path.append(source_location)
from scripts.utils import get_sabdab_details

from scripts.abag_interactions_hydrophobic import *
from scripts.abag_interactions_rings import *
from scripts.more_utils import *
AA_LIST = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
           "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
casa_dir = source_location
data_dir = Path.joinpath(casa_dir, "data")
str_dir = Path.joinpath(casa_dir, "structures/raw")
exposed_dir = Path.joinpath(casa_dir, "structures/exposed")
# Many PDBs have no interacting light chains, so statistics realitve to the
# proportion of their interacting polar atoms are skewed since many polar atoms
# show as noninteracting. 
def count_ONs_res(topologia, resSeq, chainID):
    cnt_SC = 0
    cnt_BB = 0
    
    for chain in topologia.chains:
        if chain.chain_id == chainID:
            for res in chain.residues:
                if res.resSeq == resSeq:
                    for atm in res.atoms:
                        elemento = atm.element.symbol
                        is_polar = elemento == 'O' or elemento == 'N'
                        if is_polar and atm.is_sidechain:
                            cnt_SC+= 1
                        elif is_polar and (not atm.is_sidechain):
                            cnt_BB+= 1
    return cnt_SC, cnt_BB

if __name__ == '__main__':

    with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
        filenames = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
        chains = pickle.load(file)
    with open(Path.joinpath(
        casa_dir, "data", 'buried_ab_ag_interface_res.pickle'), 'rb') as file:
        df_interface = pickle.load(file)

    pdb_list = list(filenames.keys())
    df_dataset = get_df_dataset(casa_dir)
    count_polar = {}
    for pdb_idcode in pdb_list:
        print(f"{pdb_idcode}", flush=True)
        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
        
        ###
        # ANTIBODY
        ###
        interface_ab = []
        for linea in df_interface.query(f"idcode == '{pdb_idcode}'").ab_ag_interface_res.values[0]:
            chainID = linea[0]
            resSeq = linea[1]
            resname = linea[2]
            # Including resname in the mix due to heavychain CDR3 extra residues with
            # a letter at the end. They will show up with the same resSeq in mdtraj
            # topology, but as long as I'm querying them, there's no problem.
            interface_ab.append(chainID+resname+str(resSeq))

        cnt_SC_ab = 0
        cnt_BB_ab = 0    
        for unique_ in set(interface_ab):
            chainID = unique_[0]
            resname = unique_[1:4]
            resSeq = int(unique_[4:])
            cnt_SC, cnt_BB = count_ONs_res(trj_in.topology, resSeq, chainID)
            cnt_SC_ab += cnt_SC
            cnt_BB_ab += cnt_BB

        ###
        # ANTIGEN
        ###
        interface_ag = []
        for linea in df_interface.query(f"idcode == '{pdb_idcode}'").ag_ab_interface_res.values[0]:
            chainID = linea[0]
            resSeq = linea[1]
            resname = linea[2]
            # Including resname in the mix due to heavychain CDR3 extra residues with
            # a letter at the end. They will show up with the same resSeq in mdtraj
            # topology, but as long as I'm querying them, there's no problem.
            interface_ag.append(chainID+resname+str(resSeq))

        cnt_SC_ag = 0
        cnt_BB_ag = 0
        for unique_ in set(interface_ag):
            chainID = unique_[0]
            resname = unique_[1:4]
            resSeq = int(unique_[4:])
            cnt_SC, cnt_BB = count_ONs_res(trj_in.topology, resSeq, chainID)
            cnt_SC_ag += cnt_SC
            cnt_BB_ag += cnt_BB
        
        count_polar[pdb_idcode] = PolarCount(
            cdr_SC=cnt_SC_ab, cdr_BB=cnt_BB_ab, epi_SC=cnt_SC_ag, epi_BB=cnt_BB_ag)
        
    with open(Path.joinpath(casa_dir, 'data', 'count_polar.pkl'), 'wb') as file:
        pickle.dump(count_polar, file)


  