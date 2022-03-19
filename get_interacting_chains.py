import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import Counter
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
HBondAtom = namedtuple('HBondAtom', ['chainID', 'chain_type',
                       'CDR', 'resSeq', 'resname', 'index', 'serial', 'element', 'is_sidechain'])
HBond = namedtuple('HBond', ['donor', 'acceptor'])
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
    with open(Path.joinpath(casa_dir, 'data', 'hbonds_0_39.pkl'), 'rb') as file:
        hbonds = pickle.load(file)

    pdb_list = list(filenames.keys())
    df_dataset = get_df_dataset(casa_dir)

    interacting_chains = {}
    for pdb_idcode in pdb_list:
        print(f"{pdb_idcode}", flush=True)
        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))

        serial_to_id = {}
        for atomo in trj_in.topology.atoms:
            serial_to_id[atomo.serial] = atomo.index

        ids_oxy_nitro = []
        for lista_hbond in hbonds[pdb_idcode].values():
            for hb in lista_hbond:
                ids_oxy_nitro.append(hb.donor.chainID)
                ids_oxy_nitro.append(hb.acceptor.chainID)

        interacting_chains[pdb_idcode] = set(Counter(ids_oxy_nitro).keys())

    with open(Path.joinpath(data_dir, 'interacting_chains.pkl'), 'wb') as file:
        pickle.dump(interacting_chains, file)
