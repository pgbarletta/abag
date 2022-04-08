import sys
from pathlib import Path
import pickle
import pandas as pd
import mdtraj as md
import networkx as nx
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
PiPiPair = namedtuple('PiPiPair', ['antibody', 'antigen'])
PionPair = namedtuple('PionPair', ['ring', 'ion'])

source_location = Path().resolve()
sys.path.append(source_location)
from scripts.utils import get_sabdab_details

from scripts.abag_interactions_hydrophobic import *
from scripts.abag_interactions_rings import *
from scripts.more_utils import *

AA_LIST = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
           "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
casa_dir = source_location
str_dir = Path.joinpath(casa_dir, "structures/raw")
exposed_dir = Path.joinpath(casa_dir, "structures/exposed")

with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
    filenames = pickle.load(file)
with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
    chains = pickle.load(file)

pdb_list = list(filenames.keys())
df_dataset = get_df_dataset(casa_dir)

############
# RUN-
############

df_hydrophobic_atom_indices = pd.DataFrame(
    columns=['idcode', 'atom_indices'])
df_hydrophobic_atom_serials = pd.DataFrame(
    columns=['idcode', 'atom_serials'])
df_hydrophobic_resSeq = pd.DataFrame(
    columns=['idcode', 'resSeq'])
df_hydrophobic_resname = pd.DataFrame(
    columns=['idcode', 'resname'])
df_hydrophobic_chain_ID = pd.DataFrame(
    columns=['idcode', 'chain_ID'])
df_hydrophobic_chain_type = pd.DataFrame(
    columns=['idcode', 'chain_type'])
df_hydrophobic_cdr = pd.DataFrame(
    columns=['idcode', 'CDR'])

# check_pdb = '1qfw'
# idx = pdb_list.index(check_pdb)
# for pdb_idcode in [pdb_list[idx]]:
for pdb_idcode in pdb_list:
    print(f"{pdb_idcode}", flush=True)

    pdb_filename = Path(filenames[pdb_idcode])
    trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
    ab_chains = chains[pdb_idcode].antibody
    ag_chain = chains[pdb_idcode].antigen
    #
    # Hydrophobic clusters
    #
    try:
        G = get_carbons_graph(trj_in, df_dataset, pdb_idcode,
                              ab_chains, ag_chain, cutoff_carbons)
        pre_clusteres = get_putative_clusters(G)
        clusters = merge_clusters(trj_in, pre_clusteres, cutoff_clusters)
        all_atom_indices, all_atom_serials, all_resSeq, all_resname,\
            all_chainIDs, all_chain_type, all_cdr = get_data_from_clusteres(
                pdb_idcode, trj_in, clusters, ab_chains, df_dataset)
    except Exception as e:
        all_atom_indices = []
        all_atom_serials = []
        all_resSeq = []
        all_resname = []
        all_chainIDs = []
        all_chain_type = []
        all_cdr = []
        logging.warning(
            f"- {pdb_idcode} raised: {e.__class__}, saying: {e}, during hydrophobic "
            f"interactions calculation. Probably has no hydrophobic interactions.")
        # raise e
    finally:
        # Collect
        df_hydrophobic_atom_indices = pd.concat([df_hydrophobic_atom_indices, pd.DataFrame({
            'idcode': pdb_idcode, 'atom_indices': [np.array(all_atom_indices, dtype=object)]})])
        df_hydrophobic_atom_serials = pd.concat([df_hydrophobic_atom_serials, pd.DataFrame({
            'idcode': pdb_idcode, 'atom_serials': [np.array(all_atom_serials, dtype=object)]})])
        df_hydrophobic_resSeq = pd.concat([df_hydrophobic_resSeq, pd.DataFrame({
            'idcode': pdb_idcode, 'resSeq': [np.array(all_resSeq, dtype=object)]})])
        df_hydrophobic_resname = pd.concat([df_hydrophobic_resname, pd.DataFrame({
            'idcode': pdb_idcode, 'resname': [np.array(all_resname, dtype=object)]})])
        df_hydrophobic_chain_ID = pd.concat([df_hydrophobic_chain_ID, pd.DataFrame({
            'idcode': pdb_idcode, 'chain_ID': [np.array(all_chainIDs, dtype=object)]})])
        df_hydrophobic_chain_type = pd.concat([df_hydrophobic_chain_type, pd.DataFrame({
            'idcode': pdb_idcode, 'chain_type': [np.array(all_chain_type, dtype=object)]})])
        df_hydrophobic_cdr = pd.concat([df_hydrophobic_cdr, pd.DataFrame({
            'idcode': pdb_idcode, 'CDR': [np.array(all_cdr, dtype=object)]})])


lista_df_hydrophobic = (
    [df_hydrophobic_atom_indices, df_hydrophobic_atom_serials,
     df_hydrophobic_resSeq, df_hydrophobic_resname, df_hydrophobic_chain_ID,
     df_hydrophobic_chain_type, df_hydrophobic_cdr])

with open(Path.joinpath(casa_dir, 'data', '4_hydrophobic.pkl'), 'wb') as file:
    pickle.dump(lista_df_hydrophobic, file)

print(f" --- Done -- ")
