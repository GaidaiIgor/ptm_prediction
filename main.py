import sys
import os
import os.path as path
import itertools
import time
import math

import Bio.PDB as PDB
import scipy.spatial as spatial
import numpy as numpy


class AtomPoint():
    pass


class SecondaryStructureInfo():
    def __init__(self, chain_id, interval_start, interval_end, secondary_structure_class):
        self.chain_id = chain_id
        self.interval_start = interval_start
        self.interval_end = interval_end
        self.secondary_structure_class = secondary_structure_class

# noinspection PyTypeChecker
# http://blog.marmakoide.org/?p=1
def generate_sphere_points(n):
    golden_angle = numpy.pi * (3 - numpy.sqrt(5))
    theta = golden_angle * numpy.arange(n)
    z = numpy.linspace(1, - 1, n)
    radius = numpy.sqrt(1 - z * z)
    points = numpy.zeros((n, 3))
    points[:, 0] = radius * numpy.cos(theta)
    points[:, 1] = radius * numpy.sin(theta)
    points[:, 2] = z

    return points


def square_distance(point1, point2):
    return sum([(coord1 - coord2) ** 2 for (coord1, coord2) in zip(point1, point2)])


def make_boxes(points, box_side):
    boxes = dict()
    for point in points:
        box_index = tuple(map(lambda coord: int(coord / box_side), point.coord))
        if box_index in boxes:
            boxes[box_index].append(point)
        else:
            boxes[box_index] = [point]

    return boxes


def shift_index(shift, index):
    return tuple(shift + index for (shift, index) in zip(shift, index))


def collect_points(box_index, boxes):
    points = []
    for shift in itertools.product([-1, 0, 1], repeat = 3):
        next_index = shift_index(shift, box_index)
        if next_index in boxes:
            for point in boxes[next_index]:
                if not hasattr(point, 'deleted'):
                    points.append(point)

    return points


def remove_inner_points(protein_residues, solvent_radius, max_atom_radius):
    box_side = max_atom_radius + solvent_radius
    atom_points = [atom_point for residue in protein_residues for atom in residue if hasattr(atom, 'sas_points')
                   for atom_point in atom.sas_points]
    point_boxes = make_boxes(atom_points, box_side)
    atoms = [atom for residue in protein_residues for atom in residue if hasattr(atom, 'sas_points')]
    atom_boxes = make_boxes(atoms, box_side)

    for box_index in atom_boxes:
        box_neighbour_points = collect_points(box_index, point_boxes)
        for atom in atom_boxes[box_index]:
            for point in box_neighbour_points:
                if hasattr(point, 'deleted'):
                    continue
                if square_distance(atom.coord, point.coord) < (atom.radius + solvent_radius) ** 2 and point.parent != atom:
                    point.deleted = True

                    if point.id == 1 and point.parent.id == 'O' and point.parent.parent.id[1] == 497:
                        print(atom)
                        print(atom.parent)

    # now replace sas_points with sifted
    for atom in atoms:
        sas_points = []
        for point in atom.sas_points:
            if not hasattr(point, 'deleted'):
                sas_points.append(point)
        atom.sas_points = sas_points

    return protein_residues


def add_atom_surface(atom, sphere_points_amount, solvent_radius):
    sphere_radius = atom.radius + solvent_radius
    sphere_points = generate_sphere_points(sphere_points_amount) * sphere_radius + atom.get_coord()
    atom_points = []

    point_id = 0
    for point in sphere_points:
        next_atom_point = AtomPoint()
        next_atom_point.coord = list(point)
        next_atom_point.parent = atom
        next_atom_point.id = point_id
        point_id += 1
        atom_points.append(next_atom_point)

    atom.sas_points = atom_points
    atom.max_sphere_points = sphere_points_amount


def build_unmodified_sas(structure, sphere_points_amount, solvent_radius):
    # returns array of protein residues with added to standard residue atoms 'sas_points' field
    # which contains surface points of protein (only external points, inner points are removed)

    protein_residues = []
    for residue in structure.get_residues():
        # we don't need water
        if residue.id[0] == 'W':
            continue
        # if it's a modified residue
        # so we should build spheres only around standard atoms
        # because we want to predict modification that didn't happen yet
        elif residue.id[0].startswith('H_'):
            # skip non amino acid het atoms
            if residue.resname not in structure.mod_code_to_unmod_code:
                continue
            protein_residues.append(residue)
            for possible_atom_names in standard_residue_atoms[structure.mod_code_to_unmod_code[residue.resname]]:
                # create 1-element tuple to avoid iterating over string letters
                if not isinstance(possible_atom_names, tuple):
                    possible_atom_names = (possible_atom_names, )
                for atom_name in possible_atom_names:
                    if residue.has_id(atom_name):
                        atom = residue[atom_name]
                        add_atom_surface(atom, sphere_points_amount, solvent_radius)
                        break
                # residue lacks some necessary atom
                else:
                    print('file: ' + structure.id + ', chain: ' + residue.parent.id +
                          ', residue: ' + str(residue.id) + ' lacks some necessary atom')
        # if residue is not modified then just iterate over all atoms and check that every atom present in standard atoms
        else:
            protein_residues.append(residue)
            for atom in residue:
                # skip hydrogen atoms
                if atom.element == 'H':
                    continue
                # skip terminal oxygen
                if atom.id == 'OXT':
                    continue
                # check if atom is in standard residue atoms
                for possible_atom_names in standard_residue_atoms[residue.resname]:
                    if not isinstance(possible_atom_names, tuple):
                        possible_atom_names = (possible_atom_names, )
                    if atom.id in possible_atom_names:
                        add_atom_surface(atom, sphere_points_amount, solvent_radius)
                        break
                else:
                    print('file: ' + structure.id + ', chain: ' + residue.parent.id +
                          ', residue: ' + str(residue.id) + ', atom: ' + atom.id + ' is not in standard atoms')
    max_atom_radius = max(standard_atom_radius.values())
    return remove_inner_points(protein_residues, solvent_radius, max_atom_radius)


# calculates residue exposure relative to residue size
def residue_exposure(residue):
    total_points = 0
    possible_points = 0
    for atom in residue:
        # don't count backbone atoms
        if hasattr(atom, 'sas_points') and atom.id not in standard_backbone_atoms:
            total_points += len(atom.sas_points)
            possible_points += atom.max_sphere_points

    return total_points / possible_points


# returns one of existing atoms or None if no one exists
def get_any_atom(atom_names, residue):
    if not isinstance(atom_names, tuple):
        atom_names = (atom_names, )

    # can't collect data for not full residue
    for atom_name in atom_names:
        if atom_name in residue:
            return residue[atom_name]
    return None


def residue_depth(residue, surface_kdtree):
    furthest_atom_names = furthest_atoms[residue.canonical_name]
    furthest_atom = get_any_atom(furthest_atom_names, residue)

    if furthest_atom is None:
        return ['', '']
    nearest_surface_point_info = surface_kdtree.query(furthest_atom.coord)

    return nearest_surface_point_info


# returns distance theta and phi (spherical coordinates with center in base_atom) for atom
def positional_features(base_atom, atom):
    relative_atom_coord = atom.coord - base_atom.coord
    distance = sum(map(lambda x: x ** 2, relative_atom_coord)) ** .5
    theta = math.acos(relative_atom_coord[2] / distance)
    phi = math.atan(relative_atom_coord[1] / relative_atom_coord[0])
    return [distance, theta, phi]


# returns type, exposure and depth of residue
def amino_acid_features(residue, surface_kdtree):
    res_depth = residue_depth(residue, surface_kdtree)
    return [residue.canonical_name, residue_exposure(residue), res_depth[0], residue.secondary_structure_class]


def neighbours_features(residue, neighbours_amount, atoms, atom_kdtree, surface_kdtree):
    features = []
    query_size = neighbours_amount * max_residue_size
    furthest_atom_names = furthest_atoms[residue.canonical_name]
    furthest_atom = get_any_atom(furthest_atom_names, residue)

    if furthest_atom is None:
        return [''] * (len(positional_general_headers) + len(amino_acid_general_headers)) * neighbours_amount

    used_residues = {residue.id[1]}
    neighbour_atoms_info = atom_kdtree.query(furthest_atom.coord, query_size)
    for neighbour_number in neighbour_atoms_info[1]:
        atom = atoms[neighbour_number]
        parent_residue = atom.parent

        # skip atoms which belong to already described residues
        if parent_residue.id[1] in used_residues:
            continue
        used_residues.add(parent_residue.id[1])
        features += positional_features(furthest_atom, atom)
        features += amino_acid_features(parent_residue, surface_kdtree)

        if len(used_residues) - 1 >= neighbours_amount:
            break

    # shouldn't happen
    if len(used_residues) - 1 < neighbours_amount:
        print('Found less residues than required')

    return features


def collect_features(residue, neighbours_amount, sas):
    filename = residue.parent.parent.parent.id
    res_status = 'modified' if residue.canonical_name != residue.resname else 'unmodified'

    atom_points = [atom_point for residue in sas for atom in residue if hasattr(atom, 'sas_points')
                   for atom_point in atom.sas_points]
    point_coords = [point.coord for point in atom_points]
    atoms = [atom for residue in sas for atom in residue if hasattr(atom, 'sas_points')]
    atom_coords = [atom.coord for atom in atoms]

    surface_kdtree = spatial.KDTree(point_coords)
    atom_kdtree = spatial.KDTree(atom_coords)

    features = [filename, residue.id[1], res_status] + amino_acid_features(residue, surface_kdtree) + \
               neighbours_features(residue, neighbours_amount, atoms, atom_kdtree, surface_kdtree)
    str_features = [str(feature) for feature in features]

    return str_features


def add_unmodified_residue_names(structure):
    for residue in structure.get_residues():
        if residue.resname in structure.mod_code_to_unmod_code:
            residue.canonical_name = structure.mod_code_to_unmod_code[residue.resname]
        else:
            residue.canonical_name = residue.resname


def add_atom_radii(structure):
    for atom in structure.get_atoms():
        if atom.element in standard_atom_radius:
            atom.radius = standard_atom_radius[atom.element]


def add_residue_secondary_structure(structure, secondary_structure_info):
    no_secondary_structure_class = 14

    residues = {}
    # first, let's assign no structure to all
    for residue in structure.get_residues():
        residue.secondary_structure_class = no_secondary_structure_class
        residues[residue.parent.id + str(residue.id[1])] = residue

    # now let's change it for specified intervals
    for info in secondary_structure_info:
        for i in range(info.interval_start, info.interval_end + 1):
            residues[info.chain_id + str(i)].secondary_structure_class = info.secondary_structure_class


def add_structure_info(structure, pdb_file):
    pdb_file.seek(0)
    mod_code_to_unmod_code = {}
    mod_name_to_mod_code = {}
    secondary_structure_info = []
    beta_sheet_class_shift = 12

    for line in pdb_file:
        header = line[:6].strip()
        if header == 'MODRES':
            modified_name = line[12:15].strip()
            unmodified_name = line[24:27].strip()
            mod_code_to_unmod_code[modified_name] = unmodified_name
        elif header == 'HETNAM':
            modres_code = line[11:14].strip()
            modres_full_name = line[15:].strip()
            mod_name_to_mod_code[modres_full_name] = modres_code
        elif header == 'HELIX':
            chain_id = line[19]
            interval_start = int(line[21:25].strip())
            interval_end = int(line[33:37].strip())
            helix_class = int(line[38:40].strip())
            secondary_structure_info.append(SecondaryStructureInfo(chain_id, interval_start, interval_end, helix_class))
        elif header == 'SHEET':
            chain_id = line[21]
            interval_start = int(line[22:26].strip())
            interval_end = int(line[33:37].strip())
            beta_sheet_class = int(line[38:40].strip()) + beta_sheet_class_shift
            secondary_structure_info.append(SecondaryStructureInfo(chain_id, interval_start, interval_end, beta_sheet_class))

    structure.mod_code_to_unmod_code = mod_code_to_unmod_code
    structure.mod_name_to_mod_code = mod_name_to_mod_code

    add_unmodified_residue_names(structure)
    add_residue_secondary_structure(structure, secondary_structure_info)
    add_atom_radii(structure)


# only modifications with hetnam equal to parent directory name will be extracted
def process_directory(data_path, neighbours_amount):
    print_step = 1
    sphere_points_amount = 20
    # https://ru.wikipedia.org/wiki/%D0%92%D0%BE%D0%B4%D0%B0#/media/File:Water_molecule_dimensions.svg
    oh_bond_length = 0.96
    # https://en.wikipedia.org/wiki/Van_der_Waals_radius
    h_van_der_waals_radius = 1.2
    solvent_radius = oh_bond_length + h_van_der_waals_radius
    pdb_parser = PDB.PDBParser(PERMISSIVE = False, QUIET = True)
    filenames = os.listdir(data_path)
    pdb_entries = len([filename for filename in filenames if filename.endswith('.pdb')])

    # filenames = filenames[:1]
    # filenames = ['1E93.pdb']

    files_proceeded = 0
    for filename in filenames:
        next_file_path = path.join(data_path, filename)
        if path.isdir(next_file_path):
            yield from process_directory(next_file_path, neighbours_amount)
        elif next_file_path.endswith('.pdb'):
            with open(next_file_path) as next_file:
                structure = pdb_parser.get_structure(filename, next_file)
                add_structure_info(structure, next_file)

            modification_full_name = path.split(data_path)[1]
            modified_residue_code = structure.mod_name_to_mod_code[modification_full_name]
            unmodified_residue_code = structure.mod_code_to_unmod_code[modified_residue_code]
            sas = build_unmodified_sas(structure, sphere_points_amount, solvent_radius)

            for residue in sas:
                if residue.canonical_name == unmodified_residue_code:
                    features = collect_features(residue, neighbours_amount, sas)
                    yield features
            files_proceeded += 1
        # skip non-directory and non-pdb files
        else:
            continue
        if files_proceeded % print_step == 0:
            print('{0} out of {1} files in {2} proceeded'.format(files_proceeded, pdb_entries, data_path))


def main():
    data_path = sys.argv[1]
    neighbours_amount = 4

    neighbours_headers = []
    for i in range(neighbours_amount):
        for header in positional_general_headers + amino_acid_general_headers:
            neighbours_headers.append(header + str(i))

    with open(path.join(data_path, 'output.csv'), 'w') as output:
        output.write(','.join(['filename', 'pos', 'status'] + amino_acid_general_headers + neighbours_headers)
                     + '\n')
        for entry in process_directory(data_path, neighbours_amount):
            output.write(','.join(entry) + '\n')


standard_backbone_atoms = ['N', 'CA', 'C', ('O', 'S')]

# variable names for the same atom (or atom replacement) listed in ()
standard_residue_atoms = {'GLY': [],
                          'ALA': ['CB'],
                          'VAL': ['CB', 'CG1', 'CG2'],
                          'LEU': ['CB', 'CG', 'CD1', 'CD2'],
                          'ILE': ['CB', 'CG1', 'CG2', 'CD1'],
                          'MET': ['CB', 'CG', ('SD', 'S', 'SE'), 'CE'],
                          'PHE': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
                          'TRP': ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
                          'PRO': ['CB', 'CG', 'CD'],
                          'SER': ['CB', 'OG'],
                          'THR': ['CB', 'OG1', 'CG2'],
                          'CYS': ['CB', 'SG'],
                          'TYR': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
                          'ASN': ['CB', 'CG', 'OD1', 'ND2'],
                          'GLN': ['CB', 'CG', 'CD', 'OE1', 'NE2'],
                          'ASP': ['CB', 'CG', 'OD1', 'OD2'],
                          'GLU': ['CB', 'CG', 'CD', 'OE1', 'OE2'],
                          'LYS': ['CB', 'CG', 'CD', 'CE', 'NZ'],
                          'ARG': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                          'HIS': ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']}

for key, value in standard_residue_atoms.items():
    value.extend(standard_backbone_atoms)
max_residue_size = max(map(len, standard_residue_atoms.values()))

furthest_atoms = {'GLY': 'CA',
                  'ALA': 'CB',
                  'VAL': ('CG1', 'CG2'),
                  'LEU': ('CD1', 'CD2'),
                  'ILE': 'CD1',
                  'MET': 'CE',
                  'PHE': 'CZ',
                  'TRP': 'CH2',
                  'PRO': 'CG',
                  'SER': 'OG',
                  'THR': ('OG1', 'CG2'),
                  'CYS': 'SG',
                  'TYR': 'OH',
                  'ASN': ('OD1', 'ND2'),
                  'GLN': ('OE1', 'NE2'),
                  'ASP': ('OD1', 'OD2'),
                  'GLU': ('OE1', 'OE2'),
                  'LYS': 'NZ',
                  'ARG': ('NH1', 'NH2'),
                  'HIS': 'NE2'}
# https://en.wikipedia.org/wiki/Van_der_Waals_radius
standard_atom_radius = {'C': 1.7, 'O': 1.52, 'N': 1.55, 'S': 1.8, 'SE': 1.9}

positional_general_headers = ['distance', 'theta', 'phi']
amino_acid_general_headers = ['type', 'solvent_exposure', 'residue_depth', 'secondary_structure']


# check SNP: 3NBJ - 634 A
# check multiple models
# test()
start = time.clock()
main()
print('{0} time elapsed'.format(time.clock() - start))