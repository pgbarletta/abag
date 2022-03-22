import sys
from pathlib import Path
import itertools
import pickle
import mdtraj as md
from collections import Counter
from collections import namedtuple
ResiCount = namedtuple('ResiCount', ['antibody', 'antigen'])
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
    count_resi_pdb = {}
    count_resi = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
    count_tot_pdb = {}
    count_tot = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
    for pdb_idcode in pdb_list:
        print(f"{pdb_idcode}", flush=True)
        pdb_filename = Path(filenames[pdb_idcode])
        trj_in = md.load(Path.joinpath(exposed_dir, pdb_idcode, pdb_filename))

        count_res_ab = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
        count_res_ag = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))

        ###
        # TOTAL
        ###
        count = Counter([residue.name for residue in trj_in.topology.residues])
        count_tot_this_pdb = dict(zip(AA_LIST, itertools.repeat(0, len(AA_LIST))))
        for key, val in count.items():
            count_tot_this_pdb[key] = val
            count_tot[key] += val
        count_tot_pdb[pdb_idcode] = count_tot_this_pdb
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

        for unique_ in set(interface_ab):
            resname = unique_[1:4]
            count_res_ab[resname] += 1

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

        for unique_ in set(interface_ag):
            resname = unique_[1:4]
            count_res_ag[resname] += 1

        count_resi_pdb[pdb_idcode] = ResiCount(antibody=count_res_ab, antigen=count_res_ag)
        for key in count_resi.keys():
            count_resi[key] += count_res_ab[key]
            count_resi[key] += count_res_ag[key]
        
    with open(Path.joinpath(casa_dir, 'data', 'count_resi_pdb.pkl'), 'wb') as file:
        pickle.dump(count_resi_pdb, file)
    with open(Path.joinpath(casa_dir, 'data', 'count_resi.pkl'), 'wb') as file:
        pickle.dump(count_resi, file)
    with open(Path.joinpath(casa_dir, 'data', 'count_tot_pdb.pkl'), 'wb') as file:
        pickle.dump(count_tot_pdb, file)
    with open(Path.joinpath(casa_dir, 'data', 'count_tot.pkl'), 'wb') as file:
        pickle.dump(count_tot, file)



  