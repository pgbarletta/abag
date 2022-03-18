import pickle
import sys
from pathlib import Path
import glob
import numpy as np
from collections import namedtuple
Chains = namedtuple('Chains', ['antibody', 'antigen'])

source_location = Path().resolve()
sys.path.append(source_location)
from scripts.more_utils import *

casa_dir = source_location
data_dir = Path.joinpath(casa_dir, "data")
str_dir = Path.joinpath(casa_dir, "structures/raw")

if __name__ == '__main__':
    pdb_list = []
    with open(Path.joinpath(casa_dir, "full_pdb.list"), 'r') as file:
        for linea in file:
            pdb_list.append(linea.strip())
    df_dataset = get_df_dataset(casa_dir)
    with open(Path.joinpath(
            casa_dir, "data", 'buried_ab_ag_interface_res.pickle'), 'rb') as file:
        df_interface = pickle.load(file)

    filename_dict = {}
    chains_dict = {}
    for pdb_idcode in pdb_list:
        target_heavy_chain = df_interface.query(
            f"idcode == '{pdb_idcode}'").chainID.values[0]
        pdb_string = pdb_idcode + "_complex_??_*.pdb"
        chains_list = glob.glob(str(Path.joinpath(
            casa_dir, "structures", "exposed", pdb_idcode, pdb_string)))
        assert len(chains_list) != 0, f"BAD: {pdb_idcode}. No files. "

        for pdb_full_path in chains_list:
            pdb_filename = Path(pdb_full_path).name
            H_chain = pdb_filename[13]
            L_chain = pdb_filename[14]
            AG_chain_a = pdb_filename[16]
            AG_chain_b = pdb_filename[17]
            if H_chain == target_heavy_chain:
                c = Chains(antibody=(H_chain, L_chain), antigen=(AG_chain_a, AG_chain_b))
                chains_dict[pdb_idcode] = c
                filename_dict[pdb_idcode] = pdb_filename
                break
        else:
            raise RuntimeError(f"BAD: {pdb_idcode}. "
                               "Available files are not included in the dataset. ")

    with open(Path.joinpath(data_dir, 'filenames.pkl'), 'wb') as file:
        pickle.dump(filename_dict, file)
    with open(Path.joinpath(data_dir, 'chains.pkl'), 'wb') as file:
        pickle.dump(chains_dict, file)
