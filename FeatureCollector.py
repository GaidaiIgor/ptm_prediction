import os
import os.path as path
import math

from Structure import ProteinStructure


class FeatureCollector():
    def __init__(self, neighbours_amount = 4, show_progress = True, print_step = 1):
        self.neighbours_amount = neighbours_amount
        self.show_progress = show_progress
        self.print_step = print_step
        self.positional_general_headers = ['distance', 'theta', 'phi']
        self.amino_acid_general_headers = \
            ['type', 'solvent_exposure', 'residue_depth', 'secondary_structure', 'min_edge_distance']

    @staticmethod
    def residue_depth(residue, surface_kdtree):
        furthest_atom = ProteinStructure.get_furthest_atom(residue)
        if furthest_atom is None:
            return ['', '']

        nearest_surface_point_info = surface_kdtree.query(furthest_atom.coord)
        return nearest_surface_point_info

    @staticmethod
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

    @staticmethod
    def get_residue_type(residue):
        return residue.canonical_name

    @staticmethod
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

    @staticmethod
    def get_secondary_structure_class(residue):
        return residue.secondary_structure_class

    # returns type, exposure, depth, secondary structure and edge distance of residue
    @staticmethod
    def amino_acid_features(residue, structure):
        residue_type = FeatureCollector.get_residue_type(residue)
        residue_exposure = FeatureCollector.residue_exposure(residue)
        res_depth = FeatureCollector.residue_depth(residue, structure.get_surface_kdtree())
        min_edge_dist = FeatureCollector.min_edge_distance(residue, structure)
        secondary_structure_class = FeatureCollector.get_secondary_structure_class(residue)
        # surface_stats = get_near_surface_stats(residue, structure.get_surface_kdtree(), )
        return [residue_type, residue_exposure, res_depth[0], secondary_structure_class, min_edge_dist]

    # returns distance theta and phi (spherical coordinates with center in base_atom) for atom
    @staticmethod
    def positional_features(base_atom, atom):
        relative_atom_coord = atom.coord - base_atom.coord
        distance = sum(map(lambda x: x ** 2, relative_atom_coord)) ** .5
        theta = math.acos(relative_atom_coord[2] / distance)
        phi = math.atan(relative_atom_coord[1] / relative_atom_coord[0])
        return [distance, theta, phi]

    def neighbours_features(self, residue, structure):
        features = []
        query_size = self.neighbours_amount * ProteinStructure.max_residue_size
        furthest_atom = ProteinStructure.get_furthest_atom(residue)
        atoms = structure.get_atoms()
        if furthest_atom is None:
            return [''] * (len(self.positional_general_headers) + len(self.amino_acid_general_headers)) * self.neighbours_amount

        used_residues = {residue.id[1]}
        neighbour_atoms_info = structure.get_atom_kdtree().query(furthest_atom.coord, query_size)
        for neighbour_number in neighbour_atoms_info[1]:
            atom = atoms[neighbour_number]
            parent_residue = atom.parent
            # skip atoms which belong to already described residues
            if parent_residue.id[1] in used_residues:
                continue
            used_residues.add(parent_residue.id[1])
            features += FeatureCollector.positional_features(furthest_atom, atom)
            features += FeatureCollector.amino_acid_features(parent_residue, structure)
            if len(used_residues) - 1 >= self.neighbours_amount:
                break

        # shouldn't happen
        if len(used_residues) - 1 < self.neighbours_amount:
            print('Found less residues than required')
        return features

    def collect_features(self, residue, structure):
        features = [structure.filename, residue.id[1], residue.status] + \
                   FeatureCollector.amino_acid_features(residue, structure) + \
                   self.neighbours_features(residue, structure)
        str_features = [str(feature) for feature in features]
        return str_features

    # only modifications with hetnam equal to parent directory name will be extracted
    def get_features_from_directory(self, data_path, **protein_structure_args):
        filenames = os.listdir(data_path)
        pdb_entries = len([filename for filename in filenames if filename.endswith('.pdb')])

        files_proceeded = 0
        for filename in filenames:
            next_file_path = path.join(data_path, filename)
            if path.isdir(next_file_path):
                yield from self.get_features_from_directory(next_file_path)
            elif next_file_path.endswith('.pdb'):
                with open(next_file_path) as next_file:
                    structure = ProteinStructure(next_file, **protein_structure_args)

                modification_full_name = path.split(data_path)[1]
                modified_residue_code = structure.mod_name_to_mod_code[modification_full_name]
                unmodified_residue_code = structure.mod_code_to_unmod_code[modified_residue_code]
                structure.build_unmodified_sas()

                for residue in structure.get_residues():
                    if residue.canonical_name == unmodified_residue_code:
                        features = self.collect_features(residue, structure)
                        yield features
                files_proceeded += 1
            # skip non-directory and non-pdb files
            else:
                continue
            if files_proceeded % self.print_step == 0:
                if self.show_progress:
                    print('{0} out of {1} files in {2} proceeded'.format(files_proceeded, pdb_entries, data_path))

    def export_features(self, output_file, directory_to_process):
        neighbours_headers = []
        for i in range(self.neighbours_amount):
            for header in self.positional_general_headers + self.amino_acid_general_headers:
                neighbours_headers.append(header + str(i))

        output_file.write(','.join(['filename', 'pos', 'status'] + self.amino_acid_general_headers + neighbours_headers) + '\n')
        for entry in self.get_features_from_directory(directory_to_process):
            output_file.write(','.join(entry) + '\n')