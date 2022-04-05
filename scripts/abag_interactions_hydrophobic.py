import logging
import numpy as np
import itertools
import mdtraj as md
import networkx as nx
import string
import random

########################################
# Some useful data and parameters
ring_atoms = {
    'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']}

ring_atoms_pi_ion = {
    'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']}

cutoff_ring = .5
cutoff_angle_pipi = .85
cutoff_angle_pion = .8
cutoff_clusters = .4
cutoff_carbons = .4
########################################


def get_waters(topologia):

    ids_wat = []
    for at in topologia.atoms:
        if at.residue.name == 'HOH':
            ids_wat.append(at.index)

    return ids_wat


def is_shielded(
        positions, C_cdr_id, C_epi_id, surrounding_ON_ids, angle_cdr_ON=.85,
        angle_epi_ON=-.2):
    C_cdr_xyz = positions[C_cdr_id, :]
    C_epi_xyz = positions[C_epi_id, :]

    vec_C_C = C_epi_xyz - C_cdr_xyz
    n_vec_C_C = vec_C_C / np.linalg.norm(vec_C_C)

    for ON_id in surrounding_ON_ids:
        ON_xyz = positions[ON_id, :]
        vec_cdr_ON = ON_xyz - C_cdr_xyz
        n_vec_cdr_ON = vec_cdr_ON / np.linalg.norm(vec_cdr_ON)
        vec_epi_ON = ON_xyz - C_epi_xyz
        n_vec_epi_ON = vec_epi_ON / np.linalg.norm(vec_epi_ON)

        if np.dot(n_vec_C_C, n_vec_cdr_ON) > angle_cdr_ON and\
                np.dot(n_vec_C_C, n_vec_epi_ON) < angle_epi_ON:
            return True, ON_id

    return False, 0


def turn_serial_to_id(topologia, atom_serials, pdb_idcode):
    is_complex = False
    serial_to_id = {}
    for atomo in topologia.atoms:
        serial_to_id[atomo.serial] = atomo.index

    ids_carbons = []
    ids_oxy_nitro = []
    for id in atom_serials:
        try:
            index = serial_to_id[id]
            elemento = topologia.atom(index).element.symbol
            if elemento == 'C':
                ids_carbons.append(index)
            elif elemento == 'O' or elemento == 'N':
                ids_oxy_nitro.append(index)
        except KeyError as e:
            # logging.warning(
            #     f"turn_serial_to_id(): {pdb_idcode} either has alternate "
            #     f"positions for atom serial: {id}, or this is not in the PDB file. ")
            is_complex = True
        except:
            print(f"Error with {pdb_idcode} at atom serial {id}")

    return ids_carbons, ids_oxy_nitro, is_complex


def get_ids_CON(topologia, df_dataset, pdb_idcode, ab_chains, ag_chains):

    all_epitope_atoms = list(
        set(
            itertools.chain(
                *
                [fila.epitope_atoms for index, fila
                 in df_dataset.query(
                     f"idcode == '{pdb_idcode}' and chainID in {ab_chains}").
                 iterrows()])))

    all_cdr_atoms = np.unique(
        list(
            itertools.chain(
                *
                [fila.cdr_atoms for index,
                 fila in df_dataset.query(
                     f"idcode == '{pdb_idcode}' and chainID in {ab_chains}").
                 iterrows()])))

    ids_C_epitope_atoms, ids_ON_epitope_atoms, is_complex_epitope = turn_serial_to_id(
        topologia, all_epitope_atoms, pdb_idcode)

    ids_C_cdr_atoms, ids_ON_cdr_atoms, is_complex_cdr = turn_serial_to_id(
        topologia, all_cdr_atoms, pdb_idcode)

    return ids_C_epitope_atoms, ids_C_cdr_atoms, ids_ON_epitope_atoms,\
        ids_ON_cdr_atoms, is_complex_epitope, is_complex_cdr


def get_carbons_graph(trj_in, df_dataset, pdb_idcode, ab_chains, ag_chain, cutoff):

    ids_C_epitope_atoms, ids_C_cdr_atoms, ids_ON_epitope_atoms, ids_ON_cdr_atoms,\
        is_complex_epitope, is_complex_cdr = get_ids_CON(
            trj_in.topology, df_dataset, pdb_idcode, ab_chains, ag_chain)

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

    G = nx.Graph()
    indices_close_C_C_distancias = np.where(C_C_distancias < cutoff)
    mask_close_C_ON_distancias = C_ON_distancias < cutoff
    shielding_atoms_serial = {}

    for i, j in zip(*indices_close_C_C_distancias):
        C_cdr_id = ids_C_cdr_atoms[i]
        C_epi_id = ids_C_epitope_atoms[j]
        surrounding_ON_ids = [
            ids_ON_atoms[i]
            for i in np.where(mask_close_C_ON_distancias[i, :])[0]]

        shielded, ON_id = is_shielded(
            trj_in.xyz[0],
            C_cdr_id, C_epi_id, surrounding_ON_ids)

        C1 = trj_in.topology.atom(C_cdr_id).serial
        C2 = trj_in.topology.atom(C_epi_id).serial
        ON = trj_in.topology.atom(ON_id).serial
        if shielded:
            shielding_atoms_serial[i] = trj_in.topology.atom(ON_id).serial
        else:
            G.add_edge(C_cdr_id, C_epi_id)

    return G


def get_putative_clusters(G):
    pre_clusteres = []
    for cluster in sorted(nx.connected_components(G), key=len, reverse=True):
        pre_clusteres.append(cluster)

    return pre_clusteres


def merge_clusters(trj_in, pre_clusteres, cutoff_clusters):
    def clusters_are_close(trj_in, cluster_1, cluster_2, cutoff_clusters):

        cluster_2_chainIDs = [(i, trj_in.topology.atom(
            i).residue.chain.chain_id) for i in cluster_2]

        for carbon in cluster_1:
            chainID_A = trj_in.topology.atom(carbon).residue.chain.chain_id
            same_chain_cluster_2 = [
                i for i, chainID in cluster_2_chainIDs
                if chainID == chainID_A]

            if len(same_chain_cluster_2) == 0:
                continue
            distancias = md.compute_distances(
                trj_in,
                list(
                    itertools.product(
                        [carbon],
                        list(same_chain_cluster_2))))
            close_clusters = np.any(distancias < cutoff_clusters)
            if close_clusters:
                return True
        return False

    H = nx.Graph()
    H.add_node(0)
    for i, clu_i in enumerate(pre_clusteres):
        for j in range(i + 1, len(pre_clusteres)):
            if clusters_are_close(
                    trj_in, clu_i, pre_clusteres[j],
                    cutoff_clusters):
                H.add_edge(i, j)
            else:
                H.add_node(i)
                H.add_node(j)

    clusteres = []
    for connected_clusters in sorted(
            nx.connected_components(H),
            key=len, reverse=True):
        new_cluster = []
        for i in connected_clusters:
            new_cluster.extend(pre_clusteres[i])
        clusteres.append(new_cluster)

    # Make sure that the clusters are sorted by size
    idx = np.flip(np.argsort([len(c) for c in clusteres]))
    sorted_clusters = []
    for i in idx:
        sorted_clusters.append(clusteres[i])

    return sorted_clusters


def draw_clusters(trj_in, df_dataset, pdb_idcode,
                  pdb_filename, ab_chains, ag_chains, clusters, filename):
    with open(f"aux/{filename}", "w") as fil:
        fil.write(f'from pymol import cmd\n\n')
        fil.write(f'cmd.load("{pdb_filename}")\n')
        fil.write(f'cmd.color("salmon", "')
        for chainID in ag_chains[0:-1]:
            fil.write(f'chain {chainID} or ')
        fil.write(f'chain {ag_chains[-1]}")\n')
        fil.write(f'cmd.color("atomic", "(not elem C)")\n\n')

        # Show epitope residues as lines
        epitope_residues = list(
            set(
                itertools.chain.from_iterable(
                    fila.epitope_residues for index, fila in df_dataset.query(
                        f"idcode == '{pdb_idcode}' and chainID in {ab_chains}").iterrows())))

        for resi in epitope_residues:
            fil.write(f'cmd.show("lines", "resi {resi} and chain {ag_chains[0]}")' + '\n')

        # Show CDR residues as lines
        cdr_residues = []
        chainIDs = []

        for _, fila in df_dataset.query(
                f"idcode == '{pdb_idcode}' and chainID in {ab_chains}").iterrows():

            cdr_residues.extend(list(range(fila.cdr_begin, fila.cdr_end + 1)))
            chainIDs.extend(list(itertools.repeat(fila.chainID,
                                                  fila.cdr_end - fila.cdr_begin + 1)))

        for resi, chainID in zip(cdr_residues, chainIDs):
            fil.write(
                f'cmd.show("lines", "resi {resi} and chain {chainID}")\n')

        # CD3 has several residues numbered '100 + letter'
        cdrH3_extra_residues = [''.join(tuple)
                                for tuple in list(
                                    itertools.product(
                                        ["100"],
                                        string.ascii_uppercase[0: 12]))]

        heavy_chainIDs = list(set([fila.chainID for index, fila in df_dataset[(
            df_dataset.idcode == pdb_idcode) & (df_dataset.chain_type == 'H')].iterrows()]))

        for resi in cdrH3_extra_residues:
            fil.write(f'cmd.show("lines", "resi {resi} and (')
            for chainID in heavy_chainIDs[0:-1]:
                fil.write(f'chain {chainID} or ')
            fil.write(f'chain {heavy_chainIDs[-1]})")\n')

        ##
        #
        ##

        for n, cluster in enumerate(clusters):
            linea = f''
            cluster_id = 'cluster_' + str(n + 1)
            fil.write(f'cmd.select("id ')
            for c in cluster:
                c_serial = trj_in.topology.atom(c).serial
                linea += f'{c_serial}+'
            fil.write(linea[0:-1])
            fil.write(f'")\n')
            fil.write(f'cmd.set_name("sele", "{cluster_id}")\n')
            fil.write(f'cmd.show("spheres", "{cluster_id}")\n')
            color = "%06x" % random.randint(0, 0xFFFFFF)
            fil.write(f'cmd.color("0x{color}", "{cluster_id}")\n')


def draw_interface(
        trj_in, df_dataset, pdb_idcode, ag_chains, filename):
    with open(f"tempo/{filename}", "w") as fil:
        fil.write(f'from pymol import cmd\n\n')
        fil.write(f'cmd.load("{pdb_idcode}.pdb")\n')
        fil.write(f'cmd.color("salmon", "')
        for chainID in ag_chains[0:-1]:
            fil.write(f'chain {chainID} or ')
        fil.write(f'chain {ag_chains[-1]}")\n')
        fil.write(f'cmd.color("atomic", "(not elem C)")\n\n')

        # Show epitope residues as lines
        epitope_residues = [
            fila.epitope_residues for index,
            fila in df_dataset[df_dataset.idcode == pdb_idcode].iterrows()]

        for resi in epitope_residues:
            if len(resi) == 0:
                continue
            resi_str = ''.join([str(r) + '+' for r in resi])[0:-1]
            fil.write(f'cmd.show("lines", "resi {resi_str} and (')
            for chainID in ag_chains[0:-1]:
                fil.write(f'chain {chainID} or ')
            fil.write(f'chain {ag_chains[-1]})")\n')

        # Show CDR residues as lines
        cdr_residues = [
            list(range(fila.cdr_begin, fila.cdr_end + 1)) for index,
            fila in
            df_dataset[(df_dataset.idcode == pdb_idcode)].iterrows()]
        chainIDs = [
            fila.chainID for index,
            fila in df_dataset[df_dataset.idcode == pdb_idcode].iterrows()]

        for resi, chainID in zip(cdr_residues, chainIDs):
            resi_str = ''.join([str(r) + '+' for r in resi])[0:-1]
            fil.write(
                f'cmd.show("lines", "resi {resi_str} and chain {chainID}")\n')

        # CD3 has several residues numbered '100 + letter'
        cdrH3_extra_residues = [''.join(tuple)
                                for tuple in list(
                                    itertools.product(
                                        ["100"],
                                        string.ascii_uppercase[0: 12]))]

        heavy_chainIDs = list(set([fila.chainID for index, fila in df_dataset[(
            df_dataset.idcode == pdb_idcode) & (df_dataset.chain_type == 'H')].iterrows()]))

        for resi in cdrH3_extra_residues:
            fil.write(f'cmd.show("lines", "resi {resi} and (')
            for chainID in heavy_chainIDs[0:-1]:
                fil.write(f'chain {chainID} or ')
            fil.write(f'chain {heavy_chainIDs[-1]})")\n')


def get_data_from_clusteres(
        pdb_idcode, trj_in, clusteres, ab_chains, df_dataset):
    all_atom_indices = []
    all_atom_serials = []
    all_resSeq = []
    all_resname = []
    all_chainIDs = []
    all_chain_type = []
    all_cdr = []
    for cluster in clusteres:
        cluster_indices = np.array([], dtype=int)
        cluster_atom_serials = np.array([], dtype=int)
        cluster_resSeq = np.array([], dtype=int)
        cluster_resname = np.array([], dtype=int)
        cluster_chainIDs = np.array([])
        cluster_chain_type = np.array([], dtype=str)
        cluster_cdr = np.array([], dtype=int)
        for carbon in cluster:

            index = int(carbon)
            cluster_indices = np.append(cluster_indices, index)

            serial = trj_in.topology.atom(carbon).serial
            cluster_atom_serials = np.append(cluster_atom_serials, serial)

            resSeq = trj_in.topology.atom(carbon).residue.resSeq
            cluster_resSeq = np.append(cluster_resSeq, resSeq)

            resname = trj_in.topology.atom(carbon).residue.name
            cluster_resname = np.append(cluster_resname, resname)

            chainID = trj_in.topology.atom(carbon).residue.chain.chain_id
            cluster_chainIDs = np.append(cluster_chainIDs, chainID)

            if chainID in ab_chains:
                # Get the chain type, heavy or light
                chain_type = (df_dataset.query(f"idcode == '{pdb_idcode}'").query(
                    f"chainID == '{chainID}'").chain_type).iloc[0]

                # Get the CDR
                cdr_resSeq_ranges = [
                    (row.cdr_begin, row.cdr_end) for i,
                    row in df_dataset.query(f"idcode == '{pdb_idcode}'").query(
                        f"chainID == '{chainID}'").iterrows()]
                cdr = [i + 1 for i in range(0, 3) if resSeq in range(
                    cdr_resSeq_ranges[i][0], cdr_resSeq_ranges[i][1] + 1)][0]
            else:
                # carbon from the antigen chain
                chain_type = ""
                cdr = -1
            cluster_chain_type = np.append(cluster_chain_type, chain_type)
            cluster_cdr = np.append(cluster_cdr, cdr)

        all_atom_indices.append(cluster_indices)
        all_atom_serials.append(cluster_atom_serials)
        all_resSeq.append(cluster_resSeq)
        all_resname.append(cluster_resname)
        all_chainIDs.append(cluster_chainIDs)
        all_chain_type.append(cluster_chain_type)
        all_cdr.append(cluster_cdr)

    return all_atom_indices, all_atom_serials, all_resSeq, all_resname,\
        all_chainIDs, all_chain_type, all_cdr
