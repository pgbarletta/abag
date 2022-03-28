import sys
from pathlib import Path
import itertools
import pickle
import pandas as pd
import mdtraj as md
import logging
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])
ShieldingAtom = namedtuple(
    'ShieldingAtom', ['chainID', 'chain_type', 'CDR', 'resSeq', 'resname', 'index',
                      'serial', 'element', 'is_sidechain'])
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
cutoff = .5


def get_chain_info(df_dataset, pdb_idcode, ab_chains, chainID, resSeq):
    cdr_resSeq_ranges = [(row.cdr_begin, row.cdr_end) for i, row in df_dataset.query(
        f"idcode == '{pdb_idcode}'").query(f"chainID == '{chainID}'").iterrows()]

    if chainID in ab_chains:
        chain_type = (df_dataset.query(f"idcode == '{pdb_idcode}'").query(
            f"chainID == '{chainID}'").chain_type).iloc[0]
        try:
            cdr = [i + 1 for i in range(0, 3) if resSeq in range(
                cdr_resSeq_ranges[i][0], cdr_resSeq_ranges[i][1] + 1)][0]
        except IndexError:
            logging.warning(
                f"{pdb_idcode} -- {chainID}:{resSeq}, does not belong to a CDR.")
            cdr = 0
    else:
        chain_type = -1
        cdr = -1
    return chain_type, cdr


if __name__ == '__main__':

    with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
        filenames = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
        chains = pickle.load(file)
    pdb_list = list(filenames.keys())
    df_dataset = get_df_dataset(casa_dir)

    shielding_dict = {}
    # check_pdb = '3ze1'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        print(f"{pdb_idcode}", flush=True)

        pdb_filename = Path.joinpath(str_dir, pdb_idcode + '.pdb')
        trj_in = md.load(pdb_filename)
        ab_chains = chains[pdb_idcode].antibody
        ag_chains = chains[pdb_idcode].antigen

        ids_C_epitope_atoms, ids_C_cdr_atoms, ids_ON_epitope_atoms, ids_ON_cdr_atoms,\
            is_complex_epitope, is_complex_cdr = get_ids_CON(
                trj_in.topology, df_dataset, pdb_idcode, ab_chains, ag_chains)

        if is_complex_epitope or is_complex_cdr:
            print(f"{pdb_idcode} has epitope atoms not found on the PDB.", flush=True)

        ids_ON_atoms = ids_ON_cdr_atoms + ids_ON_epitope_atoms

        C_ON_pairs = np.array(
            list(
                itertools.product(
                    ids_C_cdr_atoms, ids_ON_atoms)))
        C_C_pairs = np.array(
            list(
                itertools.product(
                    ids_C_cdr_atoms, ids_C_epitope_atoms)))

        C_ON_distancias = md.compute_distances(
            trj_in, C_ON_pairs).reshape(
            (len(ids_C_cdr_atoms),
                len(ids_ON_atoms)))

        C_C_distancias = md.compute_distances(
            trj_in, C_C_pairs).reshape(
            (len(ids_C_cdr_atoms),
                len(ids_C_epitope_atoms)))

        indices_close_C_C_distancias = np.where(C_C_distancias < cutoff)
        mask_close_C_ON_distancias = C_ON_distancias < cutoff
        shielding_pdb = {}

        for i, j in zip(*indices_close_C_C_distancias):
            C_cdr_id = ids_C_cdr_atoms[i]
            C_epi_id = ids_C_epitope_atoms[j]
            surrounding_ON_ids = [
                ids_ON_atoms[i]
                for i in np.where(mask_close_C_ON_distancias[i, :])[0]]

            shielded, ON_id = is_shielded(
                trj_in.xyz[0],
                C_cdr_id, C_epi_id, surrounding_ON_ids)

            if shielded:
                chainID = trj_in.topology.atom(ON_id).residue.chain.chain_id
                resSeq = trj_in.topology.atom(ON_id).residue.resSeq
                resname = trj_in.topology.atom(ON_id).residue.name
                chain_type, cdr = get_chain_info(
                    df_dataset, pdb_idcode, ab_chains, chainID, resSeq)
                serial = trj_in.topology.atom(ON_id).serial
                element = trj_in.topology.atom(ON_id).element.symbol
                is_sidechain = trj_in.topology.atom(ON_id).is_sidechain

                # Compile all the info on this shielding polar atom.
                shielding_atom = ShieldingAtom(
                    chainID=chainID, chain_type=chain_type, CDR=cdr,
                    resSeq=resSeq, resname=resname, index=ON_id,
                    serial=serial, element=element, is_sidechain=is_sidechain)

                shielding_pdb[ON_id] = shielding_atom

        shielding_dict[pdb_idcode] = shielding_pdb

    with open(Path.joinpath(casa_dir, 'data', 'shielding.pkl'), 'wb') as file:
        pickle.dump(shielding_dict, file)

    print(f" --- Done -- ")
