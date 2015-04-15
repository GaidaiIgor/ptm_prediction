import os
import os.path as path
import math

import Bio.PDB as PDB

from ProteinStructure import ProteinStructure


class FeatureCollector():
    def __init__(self, neighbours_amount = 4, show_progress = True, print_step = 1, surface_radius_precision = .1):
        self.neighbours_amount = neighbours_amount
        self.show_progress = show_progress
        self.print_step = print_step
        self.surface_radius_precision = surface_radius_precision
        self.positional_general_headers = ['distance', 'theta', 'phi']
        self.amino_acid_general_headers = ['type', 'residue_sas', 'residue_ses', 'residue_depth', 'secondary_structure',
                                           'min_edge_distance', 'outer_sphere_distance', 'inner_sphere_distance', 'torsion_omega_previous',
                                           'torsion_phi', 'torsion_psi', 'torsion_omega_next']

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
    def residue_exposure(residue):
        total_ses = 0
        total_sas = 0
        for atom in residue.atoms:
            if not ProteinStructure.is_backbone_atom(atom):
                total_sas += atom.sas
                total_ses += atom.ses
        return total_sas, total_ses

    @staticmethod
    def get_secondary_structure_class(residue):
        return residue.secondary_structure_class

    @staticmethod
    def flat_list_generator(object_to_flat):
        if hasattr(object_to_flat, '__iter__') and not isinstance(object_to_flat, str):
            for item in object_to_flat:
                yield from FeatureCollector.flat_list_generator(item)
        else:
            yield object_to_flat

    @staticmethod
    def flat_list(list_to_flat):
        return list(FeatureCollector.flat_list_generator(list_to_flat))

    @staticmethod
    def get_midpoint(point1, point2):
        return [(coord1 + coord2) / 2 for (coord1, coord2) in zip(point1, point2)]

    @staticmethod
    def get_closes_point(base_point, kdtree):
        closest_point_info = kdtree.query(base_point)
        return list(kdtree.data[closest_point_info[1]])

    @staticmethod
    def square_distance(point1, point2):
        return sum([(coord1 - coord2) ** 2 for (coord1, coord2) in zip(point1, point2)])

    @staticmethod
    def find_bounding_sphere_radius(interval_start, interval_end, surface_kdtree, eps):
        base_point = interval_start
        interval_distance = FeatureCollector.square_distance(interval_start, interval_end)
        while interval_distance > eps:
            midpoint = FeatureCollector.get_midpoint(interval_start, interval_end)
            closest_point = FeatureCollector.get_closes_point(midpoint, surface_kdtree)
            if closest_point == base_point:
                interval_start = midpoint
            else:
                interval_end = midpoint
            interval_distance = FeatureCollector.square_distance(interval_start, interval_end)

        midpoint = FeatureCollector.get_midpoint(interval_start, interval_end)
        return FeatureCollector.square_distance(midpoint, base_point) ** .5

    @staticmethod
    def protein_bounding_rect_max_side(structure):
        return max([abs(sides[1] - sides[0]) for sides in structure.edge_coords])

    def get_local_surface_stats(self, surface_kdtree, base_point, normal_vector, max_distance):
        end_point = [component * max_distance for component in normal_vector]
        external_radius = FeatureCollector.find_bounding_sphere_radius(base_point, end_point, surface_kdtree, self.surface_radius_precision)
        reverse_end_point = [-coord for coord in end_point]
        internal_radius = FeatureCollector.find_bounding_sphere_radius(base_point, reverse_end_point, surface_kdtree,
                                                                       self.surface_radius_precision)
        return external_radius, internal_radius

    def get_residue_surface_stats(self, residue, structure):
        furthest_atom = ProteinStructure.get_furthest_atom(residue)
        if furthest_atom is None:
            return ['', '']

        surface_kdtree = structure.get_surface_kdtree()
        surface = structure.get_surface()
        closest_surface_point_info = surface_kdtree.query(furthest_atom.coord)
        closest_surface_point = surface[closest_surface_point_info[1]]
        max_side = FeatureCollector.protein_bounding_rect_max_side(structure)
        return self.get_local_surface_stats(surface_kdtree, closest_surface_point.coord, closest_surface_point.normal_vector, max_side)

    @staticmethod
    def calculate_torsion_omega(current_residue, next_residue):
        atom1 = current_residue['CA'].get_vector()
        atom2 = current_residue['C'].get_vector()
        atom3 = next_residue['N'].get_vector()
        atom4 = next_residue['CA'].get_vector()
        return PDB.calc_dihedral(atom1, atom2, atom3, atom4)

    @staticmethod
    def calculate_torsion_phi(previous_residue, current_residue):
        atom1 = previous_residue['C'].get_vector()
        atom2 = current_residue['N'].get_vector()
        atom3 = current_residue['CA'].get_vector()
        atom4 = current_residue['C'].get_vector()
        return PDB.calc_dihedral(atom1, atom2, atom3, atom4)

    @staticmethod
    def calculate_torsion_psi(current_residue, next_residue):
        atom1 = current_residue['N'].get_vector()
        atom2 = current_residue['CA'].get_vector()
        atom3 = current_residue['C'].get_vector()
        atom4 = next_residue['N'].get_vector()
        return PDB.calc_dihedral(atom1, atom2, atom3, atom4)

    @staticmethod
    def get_torsion_angles(residue, structure):
        residues = structure.get_residues()
        residue_index = structure.get_residue_index(residue)
        previous_residue = None
        next_residue = None
        if residue_index != 0:
            previous_residue = residues[residue_index - 1]
        if residue_index != len(residues) - 1:
            next_residue = residues[residue_index + 1]

        angles = [''] * 4
        if previous_residue is not None:
            angles[0] = FeatureCollector.calculate_torsion_omega(previous_residue, residue)
            angles[1] = FeatureCollector.calculate_torsion_phi(previous_residue, residue)
        if next_residue is not None:
            angles[2] = FeatureCollector.calculate_torsion_psi(residue, next_residue)
            angles[3] = FeatureCollector.calculate_torsion_omega(residue, next_residue)
        return angles

    def amino_acid_features(self, residue, structure):
        residue_type = FeatureCollector.get_residue_type(residue)
        residue_exposure = FeatureCollector.residue_exposure(residue)
        res_depth = FeatureCollector.residue_depth(residue, structure.get_surface_kdtree())[0]
        min_edge_dist = FeatureCollector.min_edge_distance(residue, structure)
        secondary_structure_class = FeatureCollector.get_secondary_structure_class(residue)
        surface_stats = self.get_residue_surface_stats(residue, structure)
        torsion_angles = self.get_torsion_angles(residue, structure)

        all_features = [residue_type, residue_exposure, res_depth, secondary_structure_class, min_edge_dist, surface_stats, torsion_angles]
        return FeatureCollector.flat_list(all_features)

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
        query_size = (self.neighbours_amount + 1) * ProteinStructure.max_residue_size
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
            features += self.amino_acid_features(parent_residue, structure)
            # we have collected enough information
            if len(used_residues) - 1 >= self.neighbours_amount:
                break

        assert len(used_residues) - 1 == self.neighbours_amount, "Found less residues than required in " + structure.filename
        return features

    def collect_features(self, residue, structure):
        features = [structure.filename, residue.id[1], residue.status] +\
                   self.amino_acid_features(residue, structure) +\
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
                structure.build_ses()

                # residue = structure.get_residue_by_key(ProteinStructure.make_residue_key('B', 556))
                # self.get_near_surface_stats(residue, structure)

                for residue in structure.get_residues():
                    if residue.canonical_name == unmodified_residue_code:
                        features = self.collect_features(residue, structure)
                        yield features
                files_proceeded += 1
                if files_proceeded % self.print_step == 0:
                    if self.show_progress:
                        print('{0} out of {1} files in {2} processed'.format(files_proceeded, pdb_entries, data_path))
            # skip non-directory and non-pdb files
            else:
                continue

    def export_features(self, output_file, directory_to_process):
        neighbours_headers = []
        for i in range(self.neighbours_amount):
            for header in self.positional_general_headers + self.amino_acid_general_headers:
                neighbours_headers.append(header + str(i))

        output_file.write(','.join(['filename', 'pos', 'status'] + self.amino_acid_general_headers + neighbours_headers) + '\n')
        for entry in self.get_features_from_directory(directory_to_process):
            output_file.write(','.join(entry) + '\n')

            # def split_coords(self, points):
            # xs = [point[0] for point in points]
            # ys = [point[1] for point in points]
            # zs = [point[2] for point in points]
            #     return [xs, ys, zs]

            # def plot_points(self, points):
            # figure = pyplot.figure()
            #     axes = figure.add_subplot(111, projection = '3d')
            #     [xs, ys, zs] = self.split_coords(points)
            #     axes.scatter(xs, ys, zs)
            #     pyplot.show()

            # returns type, exposure, depth, secondary structure and edge distance of residue