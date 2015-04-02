import sys
import os
import os.path as path
import time
import math

import matplotlib.pyplot as pyplot


# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import Axes3D
import scipy.spatial as spatial

from Structure import ProteinStructure


# returns type, exposure and depth of residue
def amino_acid_features(residue, surface_kdtree, structure):
    res_depth = residue_depth(residue, surface_kdtree)
    min_edge_dist = min_edge_distance(residue, structure)
    return [residue.canonical_name, residue_exposure(residue), res_depth[0], residue.secondary_structure_class, min_edge_dist]


def collect_features(residue, neighbours_amount, structure):
    point_coords = [point.coord for point in structure.get_atom_points()]
    plot_points(point_coords)
    atom_coords = [atom.coord for atom in structure.get_atoms()]

    surface_kdtree = spatial.KDTree(point_coords)
    atom_kdtree = spatial.KDTree(atom_coords)

    features = [structure.filename, residue.id[1], residue.status] + amino_acid_features(residue, surface_kdtree, structure) + \
               neighbours_features(residue, neighbours_amount, atom_kdtree, surface_kdtree, structure)
    str_features = [str(feature) for feature in features]
    return str_features


def get_surface_around(point, surface_kdtree, minimum_points):
    surface_around = []


def min_edge_distance(residue, structure):
    furthest_atom = ProteinStructure.get_furthest_atom(residue)
    if furthest_atom is None:
        return ''

    min_distance = float('inf')
    for i in range(len(furthest_atom.coord)):
        for j in range(len(structure.edge_coords[i])):
            distance = abs(furthest_atom.coord[i] - structure.edge_coords[i][j])
            if distance < min_distance:
                min_distance = distance
    return min_distance


def neighbours_features(residue, neighbours_amount, atom_kdtree, surface_kdtree, structure):
    features = []
    query_size = neighbours_amount * ProteinStructure.max_residue_size
    furthest_atom = ProteinStructure.get_furthest_atom(residue)
    atoms = structure.get_atoms()
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
        features += amino_acid_features(parent_residue, surface_kdtree, structure)
        if len(used_residues) - 1 >= neighbours_amount:
            break

    # shouldn't happen
    if len(used_residues) - 1 < neighbours_amount:
        print('Found less residues than required')
    return features


def plot_points(points):
    figure = pyplot.figure()
    axes = figure.add_subplot(111, projection = '3d')
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    zs = [point[2] for point in points]
    axes.scatter(xs, ys, zs)
    pyplot.show()


# returns distance theta and phi (spherical coordinates with center in base_atom) for atom
def positional_features(base_atom, atom):
    relative_atom_coord = atom.coord - base_atom.coord
    distance = sum(map(lambda x: x ** 2, relative_atom_coord)) ** .5
    theta = math.acos(relative_atom_coord[2] / distance)
    phi = math.atan(relative_atom_coord[1] / relative_atom_coord[0])
    return [distance, theta, phi]


# only modifications with hetnam equal to parent directory name will be extracted
def process_directory(data_path, neighbours_amount):
    print_step = 1
    sphere_points_amount = 20
    # https://ru.wikipedia.org/wiki/%D0%92%D0%BE%D0%B4%D0%B0#/media/File:Water_molecule_dimensions.svg
    oh_bond_length = 0.96
    # https://en.wikipedia.org/wiki/Van_der_Waals_radius
    h_van_der_waals_radius = 1.2
    solvent_radius = oh_bond_length + h_van_der_waals_radius
    filenames = os.listdir(data_path)
    pdb_entries = len([filename for filename in filenames if filename.endswith('.pdb')])

    files_proceeded = 0
    for filename in filenames:
        next_file_path = path.join(data_path, filename)
        if path.isdir(next_file_path):
            yield from process_directory(next_file_path, neighbours_amount)
        elif next_file_path.endswith('.pdb'):
            with open(next_file_path) as next_file:
                structure = ProteinStructure(next_file)

            modification_full_name = path.split(data_path)[1]
            modified_residue_code = structure.mod_name_to_mod_code[modification_full_name]
            unmodified_residue_code = structure.mod_code_to_unmod_code[modified_residue_code]
            structure.build_unmodified_sas(sphere_points_amount, solvent_radius)

            for residue in structure.get_residues():
                if residue.canonical_name == unmodified_residue_code:
                    features = collect_features(residue, neighbours_amount, structure)
                    yield features
            files_proceeded += 1
        # skip non-directory and non-pdb files
        else:
            continue
        if files_proceeded % print_step == 0:
            print('{0} out of {1} files in {2} proceeded'.format(files_proceeded, pdb_entries, data_path))


def residue_depth(residue, surface_kdtree):
    furthest_atom = ProteinStructure.get_furthest_atom(residue)
    if furthest_atom is None:
        return ['', '']

    nearest_surface_point_info = surface_kdtree.query(furthest_atom.coord)
    return nearest_surface_point_info


# calculates residue exposure relative to residue size
def residue_exposure(residue):
    total_points = 0
    possible_points = 0
    for atom in residue.atoms:
        # don't count backbone atoms
        if not ProteinStructure.is_backbone_atom(atom):
            total_points += len(atom.sas_points)
            possible_points += atom.max_sphere_points
    return 0 if possible_points == 0 else total_points / possible_points


def square_distance(point1, point2):
    return sum([(coord1 - coord2) ** 2 for (coord1, coord2) in zip(point1, point2)])


def test_residue_sas_points(residue, structure, solvent_radius):
    for atom in residue.atoms:
        for point in atom.sas_points:
            for other_atom in structure.get_atoms():
                if other_atom != atom:
                    if square_distance(point.coord, other_atom.coord) < (other_atom.radius + solvent_radius) ** 2:
                        return False
    return True


def main():
    data_path = sys.argv[1]
    neighbours_amount = 4

    neighbours_headers = []
    for i in range(neighbours_amount):
        for header in positional_general_headers + amino_acid_general_headers:
            neighbours_headers.append(header + str(i))

    with open(path.join(data_path, 'output1.csv'), 'w') as output:
        output.write(','.join(['filename', 'pos', 'status'] + amino_acid_general_headers + neighbours_headers) + '\n')
        for entry in process_directory(data_path, neighbours_amount):
            output.write(','.join(entry) + '\n')

positional_general_headers = ['distance', 'theta', 'phi']
amino_acid_general_headers = ['type', 'solvent_exposure', 'residue_depth', 'secondary_structure', 'min_edge_distance']

# check SNP: 3NBJ - 634 A
# check multiple models
start = time.clock()
main()
# cProfile.run('main()', '5')
print('{0} time elapsed'.format(time.clock() - start))