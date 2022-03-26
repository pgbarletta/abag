import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
PiPiPair = namedtuple('PiPiPair', ['antibody', 'antigen'])

source_location = Path().resolve()
sys.path.append(source_location)
from scripts.utils import get_sabdab_details

from scripts.abag_interactions_hydrophobic import *
from scripts.abag_interactions_rings import *
from scripts.more_utils import *

casa_dir = Path("/home/pbarletta/labo/22/AbAgInterface")
str_dir = Path.joinpath(casa_dir, "structures/raw")
exposed_dir = Path.joinpath(casa_dir, "structures/exposed")

with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
    filenames = pickle.load(file)
with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
    chains = pickle.load(file)

pdb_list = list(filenames.keys())
df_dataset = get_df_dataset(casa_dir)

# df_dataset, simple_pdb_list, simple_file_pdb_list,\
#     complex_pdb_list, complex_file_pdb_list = get_df_dataset(casa_dir)

# is_cpx_pdb = list(itertools.repeat(False, len(simple_pdb_list))) +\
#     list(itertools.repeat(True, len(complex_pdb_list)))
# pdb_list = simple_pdb_list + complex_pdb_list
# file_pdb_list = simple_file_pdb_list + complex_file_pdb_list

############
# RUN
############

df_PiPi_atom_indices = pd.DataFrame(
    columns=['idcode', 'atom_indices'])
df_PiPi_atom_serials = pd.DataFrame(
    columns=['idcode', 'atom_serials'])
df_PiPi_resSeq = pd.DataFrame(
    columns=['idcode', 'resSeq'])
df_PiPi_resname = pd.DataFrame(
    columns=['idcode', 'resname'])
df_PiPi_chain_ID = pd.DataFrame(
    columns=['idcode', 'chain_ID'])
df_PiPi_chain_type = pd.DataFrame(
    columns=['idcode', 'chain_type'])
df_PiPi_cdr = pd.DataFrame(
    columns=['idcode', 'CDR'])

# check_pdb = '2zuq'
# idx = pdb_list.index(check_pdb)
# for pdb_idcode in [pdb_list[idx]]:
for pdb_idcode in pdb_list:
    print(f"{pdb_idcode}", flush=True)

    pdb_filename = Path(filenames[pdb_idcode])
    trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
    ab_chains = chains[pdb_idcode].antibody
    ag_chain = chains[pdb_idcode].antigen

    #
    # Pi-Pi interactions
    #
    try:
        CG_rings, CoM_rings_xyz, normal_vectors = get_ring_data(
            trj_in, ab_chains, ring_atoms)
        pipi_ring_pairs = get_pipi_interactions(
            trj_in, CG_rings, CoM_rings_xyz, normal_vectors, cutoff_ring,
            cutoff_angle_pipi)

        all_atom_indices, all_atom_serials, all_resSeq, all_resname,\
            all_chainIDs, all_chain_type, all_cdr = get_data_from_ring_ring(
                pdb_idcode, trj_in, pipi_ring_pairs, ab_chains, df_dataset)
    except Exception as e:
        all_atom_indices = []
        all_atom_serials = []
        all_resSeq = []
        all_resname = []
        all_chainIDs = []
        all_chain_type = []
        all_cdr = []
        logging.error(
            f"- {pdb_idcode} raised: {e.__class__}, saying: {e}, during PiPi "
            f"interactions calculation. Aborting.")
        raise e
    finally:
        # Collect
        df_PiPi_atom_indices = pd.concat([df_PiPi_atom_indices, pd.DataFrame({
            'idcode': pdb_idcode, 'atom_indices': [all_atom_indices]})])
        df_PiPi_atom_serials = pd.concat([df_PiPi_atom_serials, pd.DataFrame({
            'idcode': pdb_idcode, 'atom_serials': [all_atom_serials]})])
        df_PiPi_resSeq = pd.concat([df_PiPi_resSeq, pd.DataFrame({
            'idcode': pdb_idcode, 'resSeq': [all_resSeq]})])
        df_PiPi_resname = pd.concat([df_PiPi_resname, pd.DataFrame({
            'idcode': pdb_idcode, 'resname': [all_resname]})])
        df_PiPi_chain_ID = pd.concat([df_PiPi_chain_ID, pd.DataFrame({
            'idcode': pdb_idcode, 'chain_ID': [all_chainIDs]})])
        df_PiPi_chain_type = pd.concat([df_PiPi_chain_type, pd.DataFrame({
            'idcode': pdb_idcode, 'chain_type': [all_chain_type]})])
        df_PiPi_cdr = pd.concat([df_PiPi_cdr, pd.DataFrame({
            'idcode': pdb_idcode, 'CDR': [all_cdr]})])


lista_df_PiPi = (
    [df_PiPi_atom_indices, df_PiPi_atom_serials, df_PiPi_resSeq,
     df_PiPi_resname, df_PiPi_chain_ID, df_PiPi_chain_type,
     df_PiPi_cdr])

with open(Path.joinpath(casa_dir, 'data', 'PiPi.pkl'), 'wb') as file:
    pickle.dump(lista_df_PiPi, file)

print(f" --- Done -- ")
