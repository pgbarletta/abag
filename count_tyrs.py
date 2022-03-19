import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
TYRs = namedtuple('TYRs', ['heavy', 'light', 'antigen'])
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


if __name__ == '__main__':

    with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
        filenames = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
        chains = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'interacting_chains.pkl'), 'rb') as file:
        interacting_chains = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'hbonds_0_39.pkl'), 'rb') as file:
        hbonds = pickle.load(file)

    pdb_list = list(filenames.keys())
    df_dataset = get_df_dataset(casa_dir)
    count_tyrs = {}
    for pdb_idcode in pdb_list:
        print(f"{pdb_idcode}", flush=True)
        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
        ab_chains = chains[pdb_idcode].antibody
        ag_chain = chains[pdb_idcode].antigen
    
        heavy_cnt = 0
        light_cnt = 0
        ag_cnt = 0
        # In `ab_chains`, heavy chain is first, then the light one.
        # This comes from the abag naming convention in structures/exposed PDBs
        for chain in trj_in.topology.chains:
            if chain.chain_id == ab_chains[0]:
                heavy_cnt = sum([ 1 for resi in chain.residues if resi.name == 'TYR' ])
            elif chain.chain_id == ab_chains[1]:
                light_cnt = sum([ 1 for resi in chain.residues if resi.name == 'TYR' ])
            elif chain.chain_id in ag_chain:
                ag_cnt = sum([ 1 for resi in chain.residues if resi.name == 'TYR' ])

        count_tyrs[pdb_idcode] = TYRs(
            heavy=heavy_cnt, light=light_cnt, antigen=ag_cnt)
        
    with open(Path.joinpath(casa_dir, 'data', 'count_tyrs.pkl'), 'wb') as file:
        pickle.dump(count_tyrs, file)


  