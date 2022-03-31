import sys
from pathlib import Path
import mdtraj as md
import subprocess
import logging
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
casa_dir = Path("/home/pbarletta/labo/22/AbAgInterface")
str_dir = Path.joinpath(casa_dir, "structures/raw")
hbond_dir = Path.joinpath(casa_dir, "hbonds")


def is_salt_bridge(donor, acceptor):
    donor_anion = (donor.resname == 'ASP' or donor.resname == 'GLU')\
        and donor.is_sidechain
    acceptor_anion = (acceptor.resname == 'ASP' or acceptor.resname == 'GLU')\
        and acceptor.is_sidechain
    donor_cation = (donor.resname == 'LYS' or donor.resname == 'ARG')\
        and donor.is_sidechain
    acceptor_cation = (acceptor.resname == 'LYS' or acceptor.resname == 'ARG')\
        and acceptor.is_sidechain

    return ((donor_anion and acceptor_cation) or (donor_cation and acceptor_anion))


def water_donor(linea, chains):

    donor_is_wat = linea[6:9] == 'HOH'
    resname = linea[20:23]
    chainID = linea[14]
    is_aminoacid = resname in AA_LIST
    is_in_ag = chainID in chains

    return (donor_is_wat and is_aminoacid and is_in_ag)


def water_acceptor(linea, chains):

    acceptor_is_wat = linea[20:23] == 'HOH'
    resname = linea[6:9]
    chainID = linea[0]
    is_aminoacid = resname in AA_LIST
    is_in_ag = chainID in chains

    return (acceptor_is_wat and is_aminoacid and is_in_ag)


def hbplus_lines_to_set(wat_lines, role_of_water="donor"):
    wat_dict = {}
    don_range = slice(0, 14)
    acc_range = slice(14, 28)

    if role_of_water == 'donor':
        wat_range = don_range
        mol_range = acc_range
    elif role_of_water == 'acceptor':
        wat_range = acc_range
        mol_range = don_range
    else:
        raise ValueError(f"hbplus_lines_to_set(): {role_of_water} is invalid argument "
                         f"to role_of_water")

    for linea in wat_lines:
        wat = linea[wat_range]
        if wat in wat_dict:
            wat_dict[wat].append(linea[mol_range])
        else:
            wat_dict[wat] = [linea[mol_range]]

    return wat_dict


def is_interface(linea, ab_chains, ag_chains):
    ag_donor_ab_acceptor = linea[0] in ag_chains and linea[14] in ab_chains
    ab_donor_ag_acceptor = linea[0] in ab_chains and linea[14] in ag_chains

    resname_donor = linea[6:9]
    resname_acceptor = linea[20:23]
    both_aminoacid = (resname_donor in AA_LIST) and (resname_acceptor in AA_LIST)

    return ((ag_donor_ab_acceptor or ab_donor_ag_acceptor) and both_aminoacid)


def get_hbonding_atom(topologia, chainID, resSeq, name):
    for cadena in topologia.chains:
        if cadena.chain_id == chainID:
            for resi in cadena.residues:
                if resi.resSeq == resSeq:
                    for atm in resi.atoms:
                        if atm.name == name:
                            return atm.index, atm.serial,\
                                atm.element.symbol, atm.is_sidechain
    raise ValueError('Hbonding atom was not found. Probably not an actual Hbond.')


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


def get_atom_from_string(pdb_idcode, linea, topologia, df_dataset, ab_chains):
    _chainID = linea[0]
    _resSeq = int(linea[1:5])
    _resname = linea[6:9]
    _chain_type, _cdr = get_chain_info(
        df_dataset, pdb_idcode, ab_chains, _chainID, _resSeq)
    _name = linea[10:14].strip()

    _index, _serial, _element, _is_sidechain = get_hbonding_atom(
        topologia, _chainID, _resSeq, _name)

    return HBondAtom(chainID=_chainID, resSeq=_resSeq, resname=_resname,
                     chain_type=_chain_type, CDR=_cdr, index=_index, serial=_serial,
                     element=_element, is_sidechain=_is_sidechain)


def get_interfacing_waters(wat_donor_ag_lines, wat_donor_ab_lines,
                           wat_acceptor_ag_lines, wat_acceptor_ab_lines):

    wat_donor_ag_dict = hbplus_lines_to_set(wat_donor_ag_lines, "donor")
    wat_donor_ab_dict = hbplus_lines_to_set(wat_donor_ab_lines, "donor")
    wat_acceptor_ag_dict = hbplus_lines_to_set(wat_acceptor_ag_lines, "acceptor")
    wat_acceptor_ab_dict = hbplus_lines_to_set(wat_acceptor_ab_lines, "acceptor")

    # Compile instances where water is the donor
    wat_donor = wat_donor_ag_dict | wat_donor_ab_dict
    interface_don_ag_acc_ab = set(wat_donor_ag_dict.keys()) &\
        set(wat_acceptor_ab_dict.keys())
    interface_don_ag_don_ab = set(wat_donor_ag_dict.keys()) &\
        set(wat_donor_ab_dict.keys())
    interface_don = {}
    for wat in [*interface_don_ag_acc_ab, *interface_don_ag_don_ab]:
        interface_don[wat] = wat_donor[wat]

    # Compile instances where water is the acceptor
    wat_acceptor = wat_acceptor_ab_dict | wat_acceptor_ag_dict
    interface_acc_ag_don_ab = set(wat_acceptor_ag_dict.keys())\
        & set(wat_donor_ab_dict.keys())
    interface_acc_ag_acc_ab = set(wat_acceptor_ag_dict.keys())\
        & set(wat_acceptor_ab_dict.keys())
    interface_acc = {}
    for wat in [*interface_acc_ag_don_ab, *interface_acc_ag_acc_ab]:
        interface_acc[wat] = wat_acceptor[wat]

    return interface_don, interface_acc


def parse_hb2(pdb_idcode, hb2_file, df_dataset, topologia, ab_chains, ag_chains):
    hbonds_dict = {}
    saltBridge_dict = {}

    wat_donor_ag_lines = []
    wat_acceptor_ag_lines = []
    wat_donor_ab_lines = []
    wat_acceptor_ab_lines = []
    with open(hb2_file, 'r') as archivo:

        lineas = archivo.readlines()
        for i, linea in enumerate(lineas):
            if i < 8:
                continue
            if water_donor(linea, ag_chains):
                wat_donor_ag_lines.append(linea)
            elif water_donor(linea, ab_chains):
                wat_donor_ab_lines.append(linea)
            elif water_acceptor(linea, ag_chains):
                wat_acceptor_ag_lines.append(linea)
            elif water_acceptor(linea, ab_chains):
                wat_acceptor_ab_lines.append(linea)
            elif is_interface(linea, ab_chains, ag_chains):
                try:
                    donor = get_atom_from_string(
                        pdb_idcode, linea[0:14], topologia, df_dataset, ab_chains)
                    acceptor = get_atom_from_string(
                        pdb_idcode, linea[14:28], topologia, df_dataset, ab_chains)

                    if is_salt_bridge(donor, acceptor):
                        saltBridge_dict[donor.serial] = [HBond(donor=donor, acceptor=acceptor)]
                        saltBridge_dict[acceptor.serial] = [HBond(donor=donor, acceptor=acceptor)]
                    else:
                        if donor.serial in hbonds_dict:
                            hbonds_dict[donor.serial].append(HBond(donor=donor, acceptor=acceptor))
                        else:
                            hbonds_dict[donor.serial] = [HBond(donor=donor, acceptor=acceptor)]
                        if acceptor.serial in hbonds_dict:
                            hbonds_dict[acceptor.serial].append(
                                HBond(donor=donor, acceptor=acceptor))
                        else:
                            hbonds_dict[acceptor.serial] = [HBond(donor=donor, acceptor=acceptor)]
                except ValueError:
                    logging.warning(f"Could not parse hb2 line: {i} {linea}")
                    continue

    # Now deal with the waters
    interface_don, interface_acc = get_interfacing_waters(
        wat_donor_ag_lines, wat_donor_ab_lines, wat_acceptor_ag_lines, wat_acceptor_ab_lines)

    for wat_str, res_list in interface_don.items():
        for res_str in res_list:
            donor = get_atom_from_string(
                pdb_idcode, wat_str, topologia, df_dataset, ab_chains)
            acceptor = get_atom_from_string(
                pdb_idcode, res_str, topologia, df_dataset, ab_chains)

            # Only add the protein residue as key
            if acceptor.serial in hbonds_dict:
                hbonds_dict[acceptor.serial].append(HBond(donor=donor, acceptor=acceptor))
            else:
                hbonds_dict[acceptor.serial] = [HBond(donor=donor, acceptor=acceptor)]

    for wat_str, res_list in interface_acc.items():
        for res_str in res_list:
            donor = get_atom_from_string(
                pdb_idcode, res_str, topologia, df_dataset, ab_chains)
            acceptor = get_atom_from_string(
                pdb_idcode, wat_str, topologia, df_dataset, ab_chains)
            # Only add the protein residue as key
            if donor.serial in hbonds_dict:
                hbonds_dict[donor.serial].append(HBond(donor=donor, acceptor=acceptor))
            else:
                hbonds_dict[donor.serial] = [HBond(donor=donor, acceptor=acceptor)]

    return hbonds_dict, saltBridge_dict


if __name__ == "__main__":
    with open(Path.joinpath(casa_dir, 'data', 'filenames.pkl'), 'rb') as file:
        filenames = pickle.load(file)
    with open(Path.joinpath(casa_dir, 'data', 'chains.pkl'), 'rb') as file:
        chains = pickle.load(file)

    pdb_list = list(filenames.keys())
    df_dataset = get_df_dataset(casa_dir)
    hbonds_dict = {}
    saltBridge_dict = {}

    # check_pdb = '4wfe'
    # idx = pdb_list.index(check_pdb)
    # for pdb_idcode in [pdb_list[idx]]:
    for pdb_idcode in pdb_list:
        print(f"{pdb_idcode}", flush=True)
        pdb_filename = Path.joinpath(str_dir, pdb_idcode + '.pdb')
        trj_in = md.load(pdb_filename)
        ab_chains = chains[pdb_idcode].antibody
        ag_chains = chains[pdb_idcode].antigen

        hbond_dir = Path.joinpath(casa_dir, "hbonds")
        hbplus = Path.joinpath(casa_dir, 'hbonds', 'hbplus')

        process = subprocess.run([hbplus, pdb_filename, "-A", "0", "0", "0", "-d", "3.9"],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, cwd=hbond_dir)

        hb2_file = Path.joinpath(hbond_dir, pdb_filename.name[0:-3] + "hb2")

        hbonds_dict[pdb_idcode], saltBridge_dict[pdb_idcode] = parse_hb2(
            pdb_idcode, hb2_file, df_dataset, trj_in.topology, ab_chains, ag_chains)

    with open(Path.joinpath(casa_dir, 'data', 'polar_hbonds.pkl'), 'wb') as file:
        pickle.dump(hbonds_dict, file)
    with open(Path.joinpath(casa_dir, 'data', 'polar_saltBridge.pkl'), 'wb') as file:
        pickle.dump(saltBridge_dict, file)
