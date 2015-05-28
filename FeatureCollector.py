import os
import os.path as path
import math

import Bio.PDB as PDB

from ProteinStructure import ProteinStructure


class FeatureCollector:
    def __init__(self, neighbours_amount=4, show_progress=True, print_step=1, surface_radius_precision=.1,
                 use_tertiary=True, protein_structure_args={}):

        self.neighbours_amount = neighbours_amount
        self.show_progress = show_progress
        self.print_step = print_step
        self.surface_radius_precision = surface_radius_precision
        self.positional_general_headers = ['distance', 'theta', 'phi']
        self.amino_acid_general_headers = ['type', 'sas', 'ses', 'depth', 'sec_struct',
                                           'edge_dist', 'outer_sphere', 'inner_sphere',
                                           'omega_prev', 'phi_tors', 'psi_tors', 'omega_next',
                                           'chain_size', 'pos', 'rel_pos']
        self.use_tertiary = use_tertiary
        self.protein_structure_args = protein_structure_args

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
    def euclidian_distance(point1, point2):
        return sum([(coord1 - coord2) ** 2 for (coord1, coord2) in zip(point1, point2)]) ** .5

    @staticmethod
    def find_bounding_sphere_radius(interval_start, interval_end, surface_kdtree, eps):
        base_point = interval_start
        interval_distance = FeatureCollector.euclidian_distance(interval_start, interval_end)
        while interval_distance > eps:
            midpoint = FeatureCollector.get_midpoint(interval_start, interval_end)
            closest_point = FeatureCollector.get_closes_point(midpoint, surface_kdtree)
            if closest_point == base_point:
                interval_start = midpoint
            else:
                interval_end = midpoint
            interval_distance = FeatureCollector.euclidian_distance(interval_start, interval_end)

        midpoint = FeatureCollector.get_midpoint(interval_start, interval_end)
        return FeatureCollector.euclidian_distance(midpoint, base_point)

    @staticmethod
    def protein_bounding_rect_max_side(structure):
        return max([abs(sides[1] - sides[0]) for sides in structure.edge_coords])

    def get_local_surface_stats(self, surface_kdtree, base_point, normal_vector, max_distance):
        end_point = [component * max_distance for component in normal_vector]
        external_radius = FeatureCollector.find_bounding_sphere_radius(base_point, end_point, surface_kdtree,
                                                                       self.surface_radius_precision)
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
        return self.get_local_surface_stats(surface_kdtree, closest_surface_point.coord, closest_surface_point.normal_vector,
                                            max_side)

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
        key_atoms = ["N", "CA", "C"]
        if residue_index != 0:
            previous_residue = residues[residue_index - 1]
            if not ProteinStructure.atoms_are_in_residue(key_atoms, previous_residue):
                previous_residue = None
        if residue_index != len(residues) - 1:
            next_residue = residues[residue_index + 1]
            if not ProteinStructure.atoms_are_in_residue(key_atoms, next_residue):
                next_residue = None

        angles = [''] * 4
        if not ProteinStructure.atoms_are_in_residue(key_atoms, residue):
            return angles
        if previous_residue is not None:
            angles[0] = FeatureCollector.calculate_torsion_omega(previous_residue, residue)
            angles[1] = FeatureCollector.calculate_torsion_phi(previous_residue, residue)
        if next_residue is not None:
            angles[2] = FeatureCollector.calculate_torsion_psi(residue, next_residue)
            angles[3] = FeatureCollector.calculate_torsion_omega(residue, next_residue)
        return angles

    @staticmethod
    def get_chain_size(residue):
        return len(residue.parent.child_list)

    @staticmethod
    def get_relative_pos(residue):
        return residue.id[1] / FeatureCollector.get_chain_size(residue)

    @staticmethod
    def get_pos(residue):
        return residue.id[1]

    def amino_acid_features(self, residue, structure):
        residue_type = FeatureCollector.get_residue_type(residue)
        residue_exposure = FeatureCollector.residue_exposure(residue)
        res_depth = FeatureCollector.residue_depth(residue, structure.get_surface_kdtree())[0]
        min_edge_dist = FeatureCollector.min_edge_distance(residue, structure)
        secondary_structure_class = FeatureCollector.get_secondary_structure_class(residue)
        surface_stats = self.get_residue_surface_stats(residue, structure)
        torsion_angles = self.get_torsion_angles(residue, structure)
        chain_size = FeatureCollector.get_chain_size(residue)
        abs_pos = FeatureCollector.get_pos(residue)
        relative_pos = FeatureCollector.get_relative_pos(residue)

        all_features = [residue_type, residue_exposure, res_depth, secondary_structure_class, min_edge_dist, surface_stats,
                        torsion_angles, chain_size, abs_pos, relative_pos]
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
        features = [structure.filename, residue.status] +\
                   self.amino_acid_features(residue, structure) +\
                   self.neighbours_features(residue, structure)
        str_features = [str(feature) for feature in features]
        return str_features

    def amino_acid_features_primary(self, residue, structure):
        residue_index = structure.get_residue_index(residue)
        residues = structure.get_residues()
        features = []
        for shift in list(range(int(-self.neighbours_amount / 2), 0)) +\
                list(range(1, int(self.neighbours_amount / 2 + 1))):
            neighbour_index = residue_index + shift
            if neighbour_index < 0 or neighbour_index >= len(residues):
                features.append("")
                continue

            neighbour = residues[neighbour_index]
            if neighbour.parent.id != residue.parent.id or neighbour.id[1] - shift != residue.id[1]:
                features.append("")
                continue
            features.append(neighbour.canonical_name)
        return features

    def collect_features_primary(self, residue, structure):
        features = [structure.filename, residue.status] + self.amino_acid_features_primary(residue, structure)
        str_features = [str(feature) for feature in features]
        return str_features

    # only modifications with hetnam equal to parent directory name will be extracted
    def get_features_from_directory(self, data_path, modification_full_name_initial):
        file_names = os.listdir(data_path)
        pdb_entries = len([filename for filename in file_names if filename.endswith('.pdb')])

        files_processed = 0
        for filename in file_names:
            next_file_path = path.join(data_path, filename)
            if path.islink(next_file_path):
                yield from self.get_features_from_directory(os.readlink(next_file_path), modification_full_name_initial)
            elif next_file_path.endswith('.pdb'):
                with open(next_file_path) as next_file:
                    structure = ProteinStructure(next_file, **self.protein_structure_args)

                if modification_full_name_initial == "auto":
                    modification_full_name = path.split(data_path)[1]
                else:
                    modification_full_name = modification_full_name_initial

                modified_residue_code = structure.mod_name_to_mod_code[modification_full_name]
                unmodified_residue_code = structure.mod_code_to_unmod_code[modified_residue_code]
                if not structure.build_ses():
                    print("Failed to build ses for {}".format(filename))
                    files_processed += 1
                    continue

                for residue in structure.get_residues():
                    if residue.canonical_name == unmodified_residue_code:
                        residue.status = "modified" if residue.resname == modified_residue_code else "unmodified"
                        if self.use_tertiary:
                            features = self.collect_features(residue, structure)
                        else:
                            features = self.collect_features_primary(residue, structure)
                        yield features
                files_processed += 1
                if files_processed % self.print_step == 0:
                    if self.show_progress:
                        print('{0} out of {1} files in {2} processed. {3}'.format(files_processed, pdb_entries,
                                                                                  data_path, filename))
            # skip non-directory and non-pdb files
            else:
                continue

    def get_header(self):
        neighbours_headers = []
        for i in range(self.neighbours_amount):
            for header in self.positional_general_headers + self.amino_acid_general_headers:
                neighbours_headers.append(header + str(i))
        return ','.join(['filename', 'status'] + self.amino_acid_general_headers + neighbours_headers) + '\n'

    def get_header_primary(self):
        neighbours_headers = []
        for i in range(self.neighbours_amount):
            neighbours_headers.append("type" + str(i))
        return ','.join(['filename', 'status'] + neighbours_headers) + '\n'

    def export_features(self, output_file, directory_to_process, modification_full_name):
        output_file.write(self.get_header() if self.use_tertiary else self.get_header_primary())
        for entry in self.get_features_from_directory(directory_to_process, modification_full_name):
            output_file.write(','.join(entry) + '\n')
