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
        pdb_string = pdb_idcode + "_complex_??_?.pdb"
        pdbs = glob.glob(str(Path.joinpath(
            casa_dir, "structures", "exposed", pdb_idcode, pdb_string)))
        if len(pdbs) == 0:
            print(f"BAD: {pdb_idcode}. No files. ", flush=True)
        conteo = np.array([sum(1 for line in open(pdb)) for pdb in pdbs])
        largest = (-conteo).argsort()
        for idx in largest:
            pdb_filename = Path(pdbs[idx]).name
            c = Chains(antibody=(pdb_filename[13], pdb_filename[14]),
                    antigen=pdb_filename[16])
            if len(df_dataset.query(f"idcode == '{pdb_idcode}' and chainID in {c.antibody}")) == 0:
                # This file is the largest, but its not included in our `df_dataset`
                continue

            chains_dict[pdb_idcode] = c
            filename_dict[pdb_idcode] = pdb_filename
            break
        else:
            print(
                f"BAD: {pdb_idcode}. Available files are not included in the dataset. ",
                flush=True)


    with open(Path.joinpath(data_dir, 'filenames.pkl'), 'wb') as file:
        pickle.dump(filename_dict, file)
    with open(Path.joinpath(data_dir, 'chains.pkl'), 'wb') as file:
        pickle.dump(chains_dict, file)
