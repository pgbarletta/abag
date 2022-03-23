import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
PionPair = namedtuple('PionPair', ['ring', 'ion'])

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

df_PiAnion_atom_indices = pd.DataFrame(
    columns=['idcode', 'atom_indices'])
df_PiAnion_atom_serials = pd.DataFrame(
    columns=['idcode', 'atom_serials'])
df_PiAnion_resSeq = pd.DataFrame(
    columns=['idcode', 'resSeq'])
df_PiAnion_resname = pd.DataFrame(
    columns=['idcode', 'resname'])
df_PiAnion_chain_ID = pd.DataFrame(
    columns=['idcode', 'chain_ID'])
df_PiAnion_chain_type = pd.DataFrame(
    columns=['idcode', 'chain_type'])
df_PiAnion_cdr = pd.DataFrame(
    columns=['idcode', 'CDR'])
df_PiCation_atom_indices = pd.DataFrame(
    columns=['idcode', 'atom_indices'])
df_PiCation_atom_serials = pd.DataFrame(
    columns=['idcode', 'atom_serials'])
df_PiCation_resSeq = pd.DataFrame(
    columns=['idcode', 'resSeq'])
df_PiCation_resname = pd.DataFrame(
    columns=['idcode', 'resname'])
df_PiCation_chain_ID = pd.DataFrame(
    columns=['idcode', 'chain_ID'])
df_PiCation_chain_type = pd.DataFrame(
    columns=['idcode', 'chain_type'])
df_PiCation_cdr = pd.DataFrame(
    columns=['idcode', 'CDR'])

# check_pdb = '3ab0'
# idx = pdb_list.index(check_pdb)
# for pdb_idcode in [pdb_list[idx]]:
for pdb_idcode in pdb_list:
    print(f"{pdb_idcode}", flush=True)

    pdb_filename = Path(filenames[pdb_idcode])
    trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))
    ab_chains = chains[pdb_idcode].antibody
    ag_chain = chains[pdb_idcode].antigen

    #
    # Pi-ion interactions
    #
    try:
        CG_rings, CoM_rings_xyz, normal_vectors = get_ring_data(
            trj_in, ab_chains, ring_atoms_pi_ion)

        _, _, ids_ON_epitope_atoms, ids_ON_cdr_atoms, _, _ = get_ids_CON(
            trj_in.topology, df_dataset, pdb_idcode, ab_chains, ag_chain)

        ring_ab_anion_ag, ring_ab_cation_ag = get_ion_ring_interactions(
            trj_in, CG_rings["antibody"], ids_ON_epitope_atoms, CoM_rings_xyz["antibody"],
            normal_vectors["antibody"], trj_in.xyz[0], cutoff_ring, cutoff_angle_pion)

        ring_ag_anion_ab, ring_ag_cation_ab = get_ion_ring_interactions(
            trj_in, CG_rings["antigen"], ids_ON_cdr_atoms, CoM_rings_xyz["antigen"],
            normal_vectors["antigen"], trj_in.xyz[0], cutoff_ring, cutoff_angle_pion)

        all_atom_indices_ring_anion,\
            all_atom_serials_ring_anion,\
            all_resSeq_ring_anion,\
            all_resname_ring_anion,\
            all_chainIDs_ring_anion,\
            all_chain_type_ring_anion,\
            all_cdr_ring_anion = get_data_from_ring_ion(
                pdb_idcode, trj_in, ring_ab_anion_ag + ring_ag_anion_ab, ab_chains, df_dataset)

        all_atom_indices_ring_cation,\
            all_atom_serials_ring_cation,\
            all_resSeq_ring_cation,\
            all_resname_ring_cation,\
            all_chainIDs_ring_cation,\
            all_chain_type_ring_cation,\
            all_cdr_ring_cation = get_data_from_ring_ion(
                pdb_idcode, trj_in, ring_ab_cation_ag + ring_ag_cation_ab, ab_chains, df_dataset)

    except Exception as e:
        all_atom_indices_ring_anion = []
        all_atom_serials_ring_anion = []
        all_resSeq_ring_anion = []
        all_resname_ring_anion = []
        all_chainIDs_ring_anion = []
        all_chain_type_ring_anion = []
        all_cdr_ring_anion = []

        all_atom_indices_ring_cation = []
        all_atom_serials_ring_cation = []
        all_resSeq_ring_cation = []
        all_resname_ring_cation = []
        all_chainIDs_ring_cation = []
        all_chain_type_ring_cation = []
        all_cdr_ring_cation = []

        logging.warning(
            f"- {pdb_idcode} raised: {e.__class__}, saying: {e}, during Pi-ion "
            f"interactions calculation. Aborting.")
        raise e
    finally:
        # Collect
        df_PiAnion_atom_indices = pd.concat([df_PiAnion_atom_indices, pd.DataFrame({
            'idcode': pdb_idcode, 'atom_indices': [all_atom_indices_ring_anion]})])
        df_PiAnion_atom_serials = pd.concat([df_PiAnion_atom_serials, pd.DataFrame({
            'idcode': pdb_idcode, 'atom_serials': [all_atom_serials_ring_anion]})])
        df_PiAnion_resSeq = pd.concat([df_PiAnion_resSeq, pd.DataFrame({
            'idcode': pdb_idcode, 'resSeq': [all_resSeq_ring_anion]})])
        df_PiAnion_resname = pd.concat([df_PiAnion_resname, pd.DataFrame({
            'idcode': pdb_idcode, 'resname': [all_resname_ring_anion]})])
        df_PiAnion_chain_ID = pd.concat([df_PiAnion_chain_ID, pd.DataFrame({
            'idcode': pdb_idcode, 'chain_ID': [all_chainIDs_ring_anion]})])
        df_PiAnion_chain_type = pd.concat([df_PiAnion_chain_type, pd.DataFrame({
            'idcode': pdb_idcode, 'chain_type': [all_chain_type_ring_anion]})])
        df_PiAnion_cdr = pd.concat([df_PiAnion_cdr, pd.DataFrame({
            'idcode': pdb_idcode, 'CDR': [all_cdr_ring_anion]})])

        df_PiCation_atom_indices = pd.concat([df_PiCation_atom_indices, pd.DataFrame(
            {'idcode': pdb_idcode, 'atom_indices': [all_atom_indices_ring_cation]})])
        df_PiCation_atom_serials = pd.concat([df_PiCation_atom_serials, pd.DataFrame(
            {'idcode': pdb_idcode, 'atom_serials': [all_atom_serials_ring_cation]})])
        df_PiCation_resSeq = pd.concat([df_PiCation_resSeq, pd.DataFrame({
            'idcode': pdb_idcode, 'resSeq': [all_resSeq_ring_cation]})])
        df_PiCation_resname = pd.concat([df_PiCation_resname, pd.DataFrame({
            'idcode': pdb_idcode, 'resname': [all_resname_ring_cation]})])
        df_PiCation_chain_ID = pd.concat([df_PiCation_chain_ID, pd.DataFrame({
            'idcode': pdb_idcode, 'chain_ID': [all_chainIDs_ring_cation]})])
        df_PiCation_chain_type = pd.concat([df_PiCation_chain_type, pd.DataFrame(
            {'idcode': pdb_idcode, 'chain_type': [all_chain_type_ring_cation]})])
        df_PiCation_cdr = pd.concat([df_PiCation_cdr, pd.DataFrame({
            'idcode': pdb_idcode, 'CDR': [all_cdr_ring_cation]})])


lista_df_PiAnion = (
    [df_PiAnion_atom_indices, df_PiAnion_atom_serials, df_PiAnion_resSeq,
     df_PiAnion_resname, df_PiAnion_chain_ID, df_PiAnion_chain_type, df_PiAnion_cdr])
with open(Path.joinpath(casa_dir, 'data', 'PiAnion.pkl'), 'wb') as file:
    pickle.dump(lista_df_PiAnion, file)

lista_df_PiCation = (
    [df_PiCation_atom_indices, df_PiCation_atom_serials, df_PiCation_resSeq,
     df_PiCation_resname, df_PiCation_chain_ID, df_PiCation_chain_type, df_PiCation_cdr])
with open(Path.joinpath(casa_dir, 'data', 'PiCation.pkl'), 'wb') as file:
    pickle.dump(lista_df_PiCation, file)

print(f" --- Done -- ")
