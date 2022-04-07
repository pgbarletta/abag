import sys
import os
from pathlib import Path
import pandas as pd
import pickle
source_location = Path().resolve()
sys.path.append(source_location)
from scripts.utils import get_sabdab_details


def lacks_epitope_atoms(buried_fullab, pdb_idcode):

    n_atoms_per_chain = [
        len(each) == 0
        for each in buried_fullab[buried_fullab.idcode
                                  == pdb_idcode].epitope_atoms]
    return all(n_atoms_per_chain)


def get_df_dataset(home_dir):
    # df_sabdab_all = pd.read_csv(Path.joinpath(
    #     source_location, 'structures/sabdab_summary_all.tsv'), sep="\t")
    df_sabdab_90 = pd.read_csv(Path.joinpath(
        source_location,
        'structures/sabdab_summary_90.tsv'),
        sep="\t")

    df_buried = pd.read_pickle(Path.joinpath(source_location,
                                             'data/epitope_buried.pickle'))

    # df_interactions = pd.read_pickle(Path.joinpath(source_location,
    #                                                'data/interactions.pickle'))

    protein_antigens = df_sabdab_90.query(
        "antigen_type == antigen_type and antigen_type.str.contains('protein')",
        engine='python').drop_duplicates()
    ab_protein_antigens = set(protein_antigens.pdb.values)
    all_saddab_proteins = set(df_sabdab_90.pdb.values)
    print(
        f"SabDab protein antigen:\n"
        f"{len(ab_protein_antigens)} proteins out of {len(all_saddab_proteins)}, "
        f"{round(len(ab_protein_antigens) / len(all_saddab_proteins) * 100, 1)}%")

    ab_both_chains = set(protein_antigens.query(
        "Hchain == Hchain and Lchain == Lchain").pdb.values)
    ab_single_H_chain = set(protein_antigens.query(
        "Hchain == Hchain").pdb.values)
    ab_single_L_chain = set(protein_antigens.query(
        "Lchain == Lchain").pdb.values)

    n_ab_no_Hchain = len(ab_protein_antigens) - len(ab_single_H_chain)
    n_ab_no_Lchain = len(ab_protein_antigens) - len(ab_single_L_chain)

    print(f"All: {len(ab_protein_antigens)}\nNo Hchain: {n_ab_no_Hchain}\nNo Lchain: {n_ab_no_Lchain}\nBoth chains: {len(ab_both_chains)}")

    buried_fullab = df_buried[df_buried.idcode.isin(ab_both_chains)]
    print(
        f"Buried surfaces of {len(set(df_buried.idcode.values))} proteins\n"
        f"with both chains: {len(set(buried_fullab.idcode.values))}"
    )

    """
    # Repeating the epitope_residues column, but as a tuple instead of a list
    buried_fullab.loc[:, 'epitope_residues_tuple'] = buried_fullab.apply(
        lambda row: tuple(row.epitope_residues), axis=1)

    # Now epitope_residues_tuple is hasheable and can be used to discard duplicates
    df_dataset = buried_fullab.drop_duplicates(
        subset=['idcode', 'chain_type', 'cdr', 'cdr_seq',
                'epitope_residues_tuple'])
    """

    # From the 867 starting PDBs with both antibody and antigen, discard the
    # 'bad' ones, that is, those without epitope atoms in the `buried_fullab`
    # DataFrame. We end up with 791 PDBs as a final data set.
    """"
    bad_pdbs = []
    with open(Path.joinpath(home_dir, "bad_pdbs.list"), 'r') as file:
        for linea in file:
            bad_pdbs.append(linea.strip())

    pre_pdb_list = list(set(buried_fullab.idcode))
    pdb_list = []
    for pdb in pre_pdb_list:
        if pdb not in bad_pdbs:
            pdb_list.append(pdb)

    with open(Path.joinpath(home_dir, 'pdb.list'), 'w') as file:
        for pdb in pdb_list:
            file.write(f"{pdb}\n")
    #
    """

    # df_dataset = buried_fullab.drop_duplicates(
    #     subset=['idcode', 'chain_type', 'cdr', 'cdr_seq'])
    # print(
    #     f" From {len(buried_fullab)} rows to {len(df_dataset)} rows",
    #     flush=True)

    simple_pdb_list = []
    with open(Path.joinpath(home_dir, "simple_pdb.list"), 'r') as file:
        for linea in file:
            simple_pdb_list.append(linea.strip())
    simple_file_pdb_list = []
    with open(Path.joinpath(home_dir, "simple_file_pdb.list"), 'r') as file:
        for linea in file:
            simple_file_pdb_list.append(linea.strip())

    complex_pdb_list = []
    with open(Path.joinpath(home_dir, "complex_pdb.list"), 'r') as file:
        for linea in file:
            complex_pdb_list.append(linea.strip())
    complex_file_pdb_list = []
    with open(Path.joinpath(home_dir, "complex_file_pdb.list"), 'r') as file:
        for linea in file:
            complex_file_pdb_list.append(linea.strip())

    bad_pdbs = []
    with open(Path.joinpath(home_dir, "bad_pdbs.list"), 'r') as file:
        for linea in file:
            bad_pdbs.append(linea.strip())

    # return buried_fullab, simple_pdb_list, simple_file_pdb_list,\
    #     complex_pdb_list, complex_file_pdb_list,
    return buried_fullab
