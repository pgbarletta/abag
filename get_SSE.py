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
res_SSE = namedtuple('res_SSE', ['index', 'resSeq', 'name', 'DSSP'])
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

    pdb_list = list(filenames.keys())
    df_dataset = get_df_dataset(casa_dir)
    
    SSE = {}
    SSE_cnt = {}
    # check_pdb = '1qfw'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        print(f"{pdb_idcode}", flush=True)
        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
        ag_chain = chains[pdb_idcode].antigen
        # Get SSE of each residue
        dssp = md.compute_dssp(trj_in)[0]
        SSE_pdb = {}
        SSE_cnt_pdb = {'H': 0, 'E': 0, 'C': 0, 'NA': 0}
        for i, res in enumerate(trj_in.topology.residues):
            if res.chain.chain_id in ag_chain:
                SSE_pdb[res.resSeq] = dssp[i]
                SSE_cnt_pdb[dssp[i]] += 1
        SSE[pdb_idcode] = SSE_pdb
        SSE_cnt[pdb_idcode] = SSE_cnt_pdb
    
    with open(Path.joinpath(casa_dir, 'data', 'SSE.pkl'), 'wb') as file:
        pickle.dump(SSE, file)
    with open(Path.joinpath(casa_dir, 'data', 'SSE_count.pkl'), 'wb') as file:
        pickle.dump(SSE_cnt, file)



  