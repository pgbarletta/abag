import numpy as np
import mdtraj as md
import itertools
import random
import logging
from collections import namedtuple
PiPiPair = namedtuple('PiPiPair', ['antibody', 'antigen'])
PionPair = namedtuple('PionPair', ['ring', 'ion'])


def is_complete(r, ring_atoms):
    min_n_atoms = len(ring_atoms[r.name]) + 4
    return not (r.n_atoms < min_n_atoms)


def get_atomic_data_from_ring_residue(positions, residue, ring_atoms):
    CD2_xyz = np.empty(3)
    CG_xyz = np.empty(3)
    CoM = np.zeros(3)
    natoms_ring = 0
    for atom in residue.atoms:
        if atom.name in ring_atoms[residue.name]:
            CoM += positions[atom.index, :]
            natoms_ring += 1
            if atom.name == 'CG':
                CG_index = atom.index
                CG_xyz = positions[atom.index, :]
            elif atom.name == 'CD2':
                CD2_xyz = positions[atom.index, :]
    CoM_xyz = CoM / natoms_ring

    v_1 = CG_xyz - CoM_xyz
    v_2 = CD2_xyz - CoM_xyz
    normal_vtor_not_norm = np.cross(v_1, v_2)

    return CG_index, CoM_xyz, normal_vtor_not_norm / np.linalg.norm(
        normal_vtor_not_norm)


def get_pipi_residues_from_chain(positions, chain, ring_atoms):
    CG_rings = []
    CoM_pipis_xzy = []
    normal_vectors = []
    for residue in chain.residues:
        if residue.name in ring_atoms:
            if is_complete(residue, ring_atoms):
                CG_index, CoM_xyz, normal_vtor = get_atomic_data_from_ring_residue(
                    positions, residue, ring_atoms)
                CG_rings.append(CG_index)
                CoM_pipis_xzy.append(CoM_xyz)
                normal_vectors.append(normal_vtor)

    return CG_rings, CoM_pipis_xzy, normal_vectors


def get_ring_data(trj_in, ab_chains, ring_atoms):
    CG_rings_ab = []
    CG_rings_ag = []
    CoM_pipis_xyz_ab = []
    CoM_pipis_xyz_ag = []
    normal_vectors_ab = []
    normal_vectors_ag = []
    for chain in trj_in.topology.chains:
        if chain.chain_id in ab_chains:
            CG_rings_ab_chain, CoM_pipis_xyz_ab_chain, normal_vectors_ab_chain =\
                get_pipi_residues_from_chain(trj_in.xyz[0], chain, ring_atoms)

            CG_rings_ab.extend(CG_rings_ab_chain)
            CoM_pipis_xyz_ab.extend(CoM_pipis_xyz_ab_chain)
            normal_vectors_ab.extend(normal_vectors_ab_chain)
        else:
            CG_rings_ag_chain, CoM_pipis_xyz_ag_chain, normal_vectors_ag_chain =\
                get_pipi_residues_from_chain(trj_in.xyz[0], chain, ring_atoms)

            CG_rings_ag.extend(CG_rings_ag_chain)
            CoM_pipis_xyz_ag.extend(CoM_pipis_xyz_ag_chain)
            normal_vectors_ag.extend(normal_vectors_ag_chain)

    return ({"antibody": CG_rings_ab, "antigen": CG_rings_ag},
            {"antibody": CoM_pipis_xyz_ab, "antigen": CoM_pipis_xyz_ag},
            {"antibody": normal_vectors_ab, "antigen": normal_vectors_ag})


def get_pipi_interactions(
        trj_in, CG_rings, CoM_rings_xyz, normal_vectors, cutoff_distance=.5,
        cutoff_angle=.9):
    pipi_ring_pairs = []
    for i, (com_i, v_i) in enumerate(
        zip(CoM_rings_xyz["antibody"],
            normal_vectors["antibody"])):
        for j, (com_j, v_j) in enumerate(
            zip(CoM_rings_xyz["antigen"],
                normal_vectors["antigen"])):
            distance = np.linalg.norm(com_i - com_j)
            angle = np.dot(v_i, v_j)
            if distance < cutoff_distance and angle > cutoff_angle:
                ring_i = trj_in.topology.atom(CG_rings["antibody"][i]).residue
                ring_j = trj_in.topology.atom(CG_rings["antigen"][j]).residue

                pipi_ring_pairs.append(
                    PiPiPair(antibody=ring_i, antigen=ring_j))
    return pipi_ring_pairs


def get_ion_ring_interactions(
    trj_in, CG_rings, ids_ON_atoms, CoM_rings_xyz, normal_vectors,
        positions, cutoff_distance, cutoff_angle=.7):

    if len(CG_rings) == 0:
        print(f"No good (whole) THR/TYR/PHE residues.")
        return [], []

    ring_ON_pairs = np.array(
        list(itertools.product(CG_rings, ids_ON_atoms)))

    ring_ON_distancias = md.compute_distances(
        trj_in, ring_ON_pairs).reshape(
        (len(CG_rings),
         len(ids_ON_atoms)))

    indices_close_ON_CG = np.where(ring_ON_distancias < cutoff_distance)
    anion_ring_pairs = []
    cation_ring_pairs = []
    for i, j in zip(*indices_close_ON_CG):
        com_xyz = CoM_rings_xyz[i]
        normal = normal_vectors[i]
        ON_xyz = positions[ids_ON_atoms[j]]
        ON_vector = ON_xyz - com_xyz
        norm_ON_vector = ON_vector / np.linalg.norm(ON_vector)

        distance = np.linalg.norm(ON_xyz - com_xyz)
        angle = np.abs(np.dot(norm_ON_vector, normal))

        if distance < cutoff_distance and angle > cutoff_angle:
            ring = trj_in.topology.atom(CG_rings[i]).residue
            ion = trj_in.topology.atom(ids_ON_atoms[j])
            if ion.element.symbol == 'O':
                anion_ring_pairs.append(PionPair(ring=ring, ion=ion))
            elif ion.element.symbol == 'N':
                cation_ring_pairs.append(PionPair(ring=ring, ion=ion))
            else:
                raise ValueError
    return anion_ring_pairs, cation_ring_pairs


def get_data_from_ring_ring(
        pdb_idcode, trj_in, mdtraj_ring_ring, ab_chains, df_dataset):
    all_atom_indices = []
    all_atom_serials = []
    all_resSeq = []
    all_resname = []
    all_chainIDs = []
    all_chain_type = []
    all_cdr = []
    for par in mdtraj_ring_ring:

        par_indices = PiPiPair(
            antibody=tuple([atom.index for atom in par.antibody.atoms]),
            antigen=tuple([atom.index for atom in par.antigen.atoms]))
        par_serials = PiPiPair(
            antibody=tuple([atom.serial for atom in par.antibody.atoms]),
            antigen=tuple([atom.serial for atom in par.antigen.atoms]))
        par_resSeqs = PiPiPair(
            antibody=par.antibody.resSeq, antigen=par.antigen.resSeq)
        par_resnames = PiPiPair(
            antibody=par.antibody.name, antigen=par.antigen.name)
        # Get chainID of the the antibody ring residue and reuse it
        ab_chainID = par.antibody.chain.chain_id
        par_chainIDs = PiPiPair(antibody=par.antibody.chain.chain_id,
                                antigen=par.antigen.chain.chain_id)
        chain_type = (df_dataset.query(f"idcode == '{pdb_idcode}'").query(
            f"chainID == '{ab_chainID}'").chain_type).iloc[0]
        par_chain_type = PiPiPair(antibody=chain_type, antigen="")

        cdr_resSeq_ranges = [
            (row.cdr_begin, row.cdr_end) for i,
            row in df_dataset.query(f"idcode == '{pdb_idcode}'").query(
                f"chainID == '{ab_chainID}'").iterrows()]
        try:
            cdr = [i + 1 for i in range(0, 3) if par_resSeqs.antibody in range(
                cdr_resSeq_ranges[i][0], cdr_resSeq_ranges[i][1] + 1)][0]
        except IndexError:
            logging.warning(
                f"{pdb_idcode} -- {par.antibody.chain.chain_id}:{par.antibody}"
                f", does not belong to a CDR.")
            cdr = 0

        par_cdr = PiPiPair(antibody=cdr, antigen=-1)

        all_atom_indices.append(par_indices)
        all_atom_serials.append(par_serials)
        all_resSeq.append(par_resSeqs)
        all_resname.append(par_resnames)
        all_chainIDs.append(par_chainIDs)
        all_chain_type.append(par_chain_type)
        all_cdr.append(par_cdr)

    return all_atom_indices, all_atom_serials, all_resSeq, all_resname,\
        all_chainIDs, all_chain_type, all_cdr


def get_data_from_ring_ion(
        pdb_idcode, trj_in, mdtraj_ring_ion, ab_chains, df_dataset):

    all_atom_indices = []
    all_atom_serials = []
    all_resSeq = []
    all_resname = []
    all_chainIDs = []
    all_chain_type = []
    all_cdr = []
    for par in mdtraj_ring_ion:

        par_indices = PionPair(
            ring=tuple([atom.index for atom in par.ring.atoms]),
            ion=par.ion.index)
        par_serials = PionPair(
            ring=tuple([atom.serial for atom in par.ring.atoms]),
            ion=par.ion.serial)
        par_resSeqs = PionPair(
            ring=par.ring.resSeq, ion=par.ion.residue.resSeq)
        par_resnames = PionPair(
            ring=par.ring.name, ion=par.ion.residue.name)
        par_chainIDs = PionPair(
            ring=par.ring.chain.chain_id, ion=par.ion.residue.chain.chain_id)
        # Figure out to which molecule has the ring and which the ion
        if par_chainIDs.ring in ab_chains:

            chain_type = (df_dataset.query(f"idcode == '{pdb_idcode}'").query(
                f"chainID == '{par_chainIDs.ring}'").chain_type).iloc[0]
            par_chain_type = PionPair(ring=chain_type, ion="")

            cdr_resSeq_ranges = [
                (row.cdr_begin, row.cdr_end) for i,
                row in df_dataset.query(f"idcode == '{pdb_idcode}'").query(
                    f"chainID == '{par_chainIDs.ring}'").iterrows()]

            try:
                cdr = [i + 1 for i in range(0, 3) if par_resSeqs.ring in range(
                    cdr_resSeq_ranges[i][0], cdr_resSeq_ranges[i][1] + 1)][0]
            except IndexError:
                logging.warning(
                    f"{pdb_idcode} -- {par.ring.chain.chain_id}:{par.ring}"
                    f", does not belong to a CDR.")
                cdr = 0
            par_cdr = PionPair(ring=cdr, ion=-1)
        else:
            chain_type = (df_dataset.query(f"idcode == '{pdb_idcode}'").query(
                f"chainID == '{par_chainIDs.ion}'").chain_type).iloc[0]
            par_chain_type = PionPair(ring="", ion=chain_type)

            cdr_resSeq_ranges = [
                (row.cdr_begin, row.cdr_end) for i,
                row in df_dataset.query(f"idcode == '{pdb_idcode}'").query(
                    f"chainID == '{par_chainIDs.ion}'").iterrows()]
            try:
                cdr = [i + 1 for i in range(0, 3) if par_resSeqs.ion in range(
                    cdr_resSeq_ranges[i][0], cdr_resSeq_ranges[i][1] + 1)][0]
            except IndexError:
                logging.warning(
                    f"{pdb_idcode} -- {par.ion.residue.chain.chain_id}:{par.ion}"
                    f", does not belong to a CDR.")
                cdr = 0
            par_cdr = PionPair(ring=-1, ion=cdr)

        all_atom_indices.append(par_indices)
        all_atom_serials.append(par_serials)
        all_resSeq.append(par_resSeqs)
        all_resname.append(par_resnames)
        all_chainIDs.append(par_chainIDs)
        all_chain_type.append(par_chain_type)
        all_cdr.append(par_cdr)

    return all_atom_indices, all_atom_serials, all_resSeq, all_resname,\
        all_chainIDs, all_chain_type, all_cdr


def draw_pi_rings(trj_in, df_dataset, pdb_idcode, ring_ion_pairs, filename):
    with open(f"tempo/{filename}", "w") as fil:
        fil.write(f'from pymol import cmd\n\n')
        fil.write(f'cmd.load("{pdb_idcode}.pdb")\n')

        for n, ring_ion in enumerate(ring_ion_pairs):
            linea = f''
            ring_ion_id = 'ring_ion_' + str(n + 1)
            fil.write(f'cmd.select("id ')
            for ri in ring_ion:
                ri_serial = trj_in.topology.atom(ri).serial
                linea += f'{ri_serial}+'
            fil.write(linea[0:-1])
            fil.write(f'")\n')
            fil.write(f'cmd.set_name("sele", "{ring_ion_id}")\n')
            fil.write(f'cmd.show("spheres", "{ring_ion_id}")\n')
            color = "%06x" % random.randint(0, 0xFFFFFF)
            fil.write(f'cmd.color("0x{color}", "{ring_ion_id}")\n')
