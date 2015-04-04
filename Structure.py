from enum import Enum
import itertools
import math

from scipy.spatial import KDTree
from Bio import PDB
import numpy


class ProteinStructure():
    def __init__(self, pdb_file = None, sphere_points_amount = 30, solvent_radius = None):
        self.atom_points = []
        self.atoms = []
        self.edge_coords = []
        self.filename = ''
        self.mod_code_to_unmod_code = {}
        self.mod_name_to_mod_code = {}
        self.protein_residues = []
        self.structure = None
        self.surface_kdtree = None
        self.atom_kdtree = None
        self.sphere_points_amount = sphere_points_amount
        if solvent_radius is None:
            self.default_solvent_radius = ProteinStructure.default_solvent_radius
        if pdb_file is not None:
            self.parse_structure(pdb_file)

    def build_unmodified_sas(self):
        for atom in self.get_atoms():
            self.__add_atom_surface(atom)
        self.__remove_inner_points()

    def get_atom_points(self):
        if not self.atom_points:
            self.atom_points = [atom_point for atom in self.get_atoms() for atom_point in atom.sas_points]
        return self.atom_points

    def get_atoms(self):
        if not self.atoms:
            self.atoms = [atom for residue in self.get_residues() for atom in residue.atoms]
        return self.atoms

    @staticmethod
    def get_furthest_atom(residue):
        furthest_atom_names = ProteinStructure.furthest_atoms[residue.canonical_name]
        furthest_atom = ProteinStructure.__get_any_atom(furthest_atom_names, residue)
        return furthest_atom

    def get_residues(self):
        # return self.protein_residues.values()
        return self.protein_residues

    @staticmethod
    def is_backbone_atom(atom):
        for atom_group in ProteinStructure.standard_backbone_atoms:
            if atom.id in atom_group:
                return True
        return False

    def parse_structure(self, pdb_file):
        pdb_parser = PDB.PDBParser(PERMISSIVE = False, QUIET = True)
        self.filename = pdb_file.name
        self.structure = pdb_parser.get_structure(pdb_file.name, pdb_file)
        self.__add_structure_info(pdb_file)

    def __add_atom_radii(self):
        for atom in self.get_atoms():
            if atom.element in ProteinStructure.standard_atom_radius:
                atom.radius = ProteinStructure.standard_atom_radius[atom.element]
            # shouldn't happen
            else:
                print('{0} at {1} in {2}: atom radius is unknown'.format(atom, atom.parent.id[1], self.filename))

    def __add_atom_surface(self, atom):
        sphere_radius = atom.radius + self.default_solvent_radius
        sphere_points = ProteinStructure.generate_sphere_points(self.sphere_points_amount) * sphere_radius + atom.coord
        atom.sas_points = []

        point_id = 0
        for point in sphere_points:
            atom.sas_points.append(AtomPoint(list(point), atom, point_id))
            point_id += 1
        atom.max_sphere_points = self.sphere_points_amount

    def __add_edge_coords(self):
        self.edge_coords = list(itertools.repeat([float('inf'), float('-inf')], 3))
        for atom in self.get_atoms():
            for i in range(len(atom.coord)):
                if atom.coord[i] < self.edge_coords[i][0]:
                    self.edge_coords[i][0] = atom.coord[i]
                if atom.coord[i] > self.edge_coords[i][1]:
                    self.edge_coords[i][1] = atom.coord[i]

    def __add_residue_secondary_structure(self, secondary_structure_info):
        # first, let's assign no structure to all
        for residue in self.get_residues():
            residue.secondary_structure_class = self.no_secondary_structure_class

        # now let's change it for specified intervals
        all_residues = self.get_residues()
        key_pos_map = dict()
        for i in range(len(all_residues)):
            key_pos_map[self.__residue_to_key(all_residues[i])] = i

        for info in secondary_structure_info:
            start_key = self.__make_residue_key(info.chain_id, info.interval_start, info.insertion_code_start)
            end_key = self.__make_residue_key(info.chain_id, info.interval_end, info.insertion_code_end)
            pos_start = key_pos_map[start_key]
            pos_end = key_pos_map[end_key]
            for i in range(pos_start, pos_end + 1):
                all_residues[i].secondary_structure_class = info.secondary_structure_class

    def __add_structure_info(self, pdb_file):
        pdb_file.seek(0)
        secondary_structure_info = []

        for line in pdb_file:
            header = line[:6].strip()
            if header == 'MODRES':
                modified_name = line[12:15].strip()
                unmodified_name = line[24:27].strip()
                self.mod_code_to_unmod_code[modified_name] = unmodified_name
            elif header == 'HETNAM':
                modres_code = line[11:14].strip()
                modres_full_name = line[15:].strip()
                self.mod_name_to_mod_code[modres_full_name] = modres_code
            elif header == 'HELIX':
                chain_id = line[19]
                interval_start = int(line[21:25].strip())
                insertion_code_start = line[25]
                interval_end = int(line[33:37].strip())
                insertion_code_end = line[37]
                helix_class = int(line[38:40].strip())
                secondary_structure_info.append(SecondaryStructureInfo(chain_id, interval_start, insertion_code_start,
                                                                       interval_end, insertion_code_end, helix_class))
            elif header == 'SHEET':
                chain_id = line[21]
                interval_start = int(line[22:26].strip())
                insertion_code_start = line[26]
                interval_end = int(line[33:37].strip())
                insertion_code_end = line[37]
                beta_sheet_class = int(line[38:40].strip()) + self.beta_sheet_class_shift
                secondary_structure_info.append(SecondaryStructureInfo(chain_id, interval_start, insertion_code_start,
                                                                       interval_end, insertion_code_end, beta_sheet_class))

        self.__init_protein_residues()
        self.__add_residue_secondary_structure(secondary_structure_info)
        self.__add_atom_radii()
        self.__add_edge_coords()

    def __add_unmodified_residue_name(self, residue):
        if residue.resname in self.mod_code_to_unmod_code:
            residue.canonical_name = self.mod_code_to_unmod_code[residue.resname]
        else:
            residue.canonical_name = residue.resname

    @staticmethod
    def __collect_points(box_index, boxes):
        for shift in itertools.product([-1, 0, 1], repeat = 3):
            next_index = ProteinStructure.__shift_index(shift, box_index)
            if next_index in boxes:
                for point in boxes[next_index]:
                    if not point.deleted:
                        yield point

    # residue.canonical_name should be set
    def __cut_residue(self, residue):
        # we should exclude all modification atoms
        residue.atoms = []
        if residue.status == ResidueStatus.modified:
            for possible_atom_names in ProteinStructure.standard_residue_atoms[residue.canonical_name]:
                for atom_name in possible_atom_names:
                    if residue.has_id(atom_name):
                        atom = residue[atom_name]
                        residue.atoms.append(atom)
                        break
                # shouldn't happen
                else:
                    print('file: ' + self.structure.id + ', chain: ' + residue.parent.id +
                          ', residue: ' + str(residue.id) + ' lacks some necessary atom')
        else:
            for atom in residue:
                # skip hydrogen atoms
                if atom.element == 'H':
                    continue
                # skip terminal oxygen
                if atom.id == 'OXT':
                    continue
                residue.atoms.append(atom)
                # check if atom is in standard residue atoms
                for possible_atom_names in ProteinStructure.standard_residue_atoms[residue.canonical_name]:
                    if atom.id in possible_atom_names:
                        break
                else:
                    print('file: ' + self.structure.id + ', chain: ' + residue.parent.id +
                          ', residue: ' + str(residue.id) + ', atom: ' + atom.id + ' is not in standard atoms')

    # http://blog.marmakoide.org/?p=1
    # noinspection PyTypeChecker
    @staticmethod
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

    # returns one of existing atoms or None if no one exists
    @staticmethod
    def __get_any_atom(atom_names, residue):
        # can't collect data for not full residue
        for atom_name in atom_names:
            if atom_name in residue:
                return residue[atom_name]
        return None

    def get_atom_kdtree(self):
        if self.atom_kdtree is None:
            points = [atom.coord for atom in self.get_atoms()]
            self.atom_kdtree = KDTree(points)
        return self.atom_kdtree

    def get_surface_kdtree(self):
        if self.surface_kdtree is None:
            points = [point.coord for point in self.get_atom_points()]
            self.surface_kdtree = KDTree(points)
        return self.surface_kdtree

    def __init_protein_residues(self):
        for residue in self.structure.get_residues():
            # we don't need water
            if residue.id[0] == 'W':
                continue
            # we don't need het atom unless it's modified residue
            if residue.id[0].startswith('H_'):
                if residue.resname not in self.mod_code_to_unmod_code:
                    # skip non amino acid het atoms
                    continue
                else:
                    residue.status = ResidueStatus.modified
            else:
                residue.status = ResidueStatus.unmodified
            self.__add_unmodified_residue_name(residue)
            self.__cut_residue(residue)
            # self.protein_residues[ProteinStructure.__residue_to_key(residue)] = residue
            self.protein_residues.append(residue)

    # it's faster than check if square distance less than square radius
    @staticmethod
    def __is_inside_sphere(test_point, sphere_center, sphere_radius):
        dz = test_point[2] - sphere_center[2]
        if -sphere_radius <= dz <= sphere_radius:
            circle_radius = sphere_radius if dz == 0 else dz * math.tan(math.acos(dz / sphere_radius))
            dy = test_point[1] - sphere_center[1]
            if -circle_radius <= dy <= circle_radius:
                x_width = circle_radius if dy == 0 else dy * math.tan(math.acos(dy / circle_radius))
                dx = test_point[0] - sphere_center[0]
                if -x_width <= dx <= x_width:
                    return True
        return False

    @staticmethod
    def __make_boxes(points, box_side):
        boxes = dict()
        for point in points:
            box_index = tuple(map(lambda coord: int(coord / box_side), point.coord))
            if box_index in boxes:
                boxes[box_index].append(point)
            else:
                boxes[box_index] = [point]
        return boxes

    @staticmethod
    def __make_list_of_lists(target_dict: dict):
        for key, value in target_dict.items():
            if not isinstance(value, list):
                target_dict[key] = [value]

    @staticmethod
    def __make_residue_key(residue_chain, residue_number, insertion_code = ' '):
        return residue_chain, residue_number, insertion_code

    def __remove_inner_points(self):
        box_side = ProteinStructure.max_atom_radius + self.default_solvent_radius
        point_boxes = ProteinStructure.__make_boxes(self.get_atom_points(), box_side)
        atom_boxes = ProteinStructure.__make_boxes(self.get_atoms(), box_side)

        for box_index in atom_boxes:
            box_neighbour_points = list(ProteinStructure.__collect_points(box_index, point_boxes))
            for atom in atom_boxes[box_index]:
                for point in box_neighbour_points:
                    if point.deleted or point.parent == atom:
                        continue
                    if ProteinStructure.__is_inside_sphere(point.coord, atom.coord, atom.radius + self.default_solvent_radius):
                        point.deleted = True

        # now replace sas_points with sifted
        for atom in self.get_atoms():
            sas_points = []
            for point in atom.sas_points:
                if not point.deleted:
                    sas_points.append(point)
            atom.sas_points = sas_points
        # invalidate atom points
        self.atom_points = []

    @staticmethod
    def __residue_to_key(residue):
        return residue.parent.id, residue.id[1], residue.id[2]

    @staticmethod
    def __shift_index(shift, index):
        return tuple([shift + index for (shift, index) in zip(shift, index)])

    beta_sheet_class_shift = 12
    no_secondary_structure_class = 14
    # solvent radius = 1.2 + 0.96 where 1.2 is oh bond length in water molecule and 0.96 is van der waals radius for hydrogen
    # https://ru.wikipedia.org/wiki/%D0%92%D0%BE%D0%B4%D0%B0#/media/File:Water_molecule_dimensions.svg
    # https://en.wikipedia.org/wiki/Van_der_Waals_radius
    default_solvent_radius = 2.16

    standard_backbone_atoms = ['N', 'CA', 'C', ['O', 'S']]
    # variable names for the same atom (or atom replacement) listed in ()
    standard_residue_atoms = {'GLY': [], 'ALA': ['CB'], 'VAL': ['CB', 'CG1', 'CG2'], 'LEU': ['CB', 'CG', 'CD1', 'CD2'],
                              'ILE': ['CB', 'CG1', 'CG2', 'CD1'], 'MET': ['CB', 'CG', ['SD', 'S', 'SE'], 'CE'],
                              'PHE': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
                              'TRP': ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
                              'PRO': ['CB', 'CG', 'CD'], 'SER': ['CB', 'OG'], 'THR': ['CB', 'OG1', 'CG2'], 'CYS': ['CB', 'SG'],
                              'TYR': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'], 'ASN': ['CB', 'CG', 'OD1', 'ND2'],
                              'GLN': ['CB', 'CG', 'CD', 'OE1', 'NE2'], 'ASP': ['CB', 'CG', 'OD1', 'OD2'],
                              'GLU': ['CB', 'CG', 'CD', 'OE1', 'OE2'], 'LYS': ['CB', 'CG', 'CD', 'CE', 'NZ'],
                              'ARG': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                              'HIS': ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2']}
    max_residue_size = max(map(len, standard_residue_atoms.values()))

    for atom_groups in standard_residue_atoms.values():
        atom_groups.extend(standard_backbone_atoms)

    furthest_atoms = {'GLY': 'CA', 'ALA': 'CB', 'VAL': ['CG1', 'CG2'], 'LEU': ['CD1', 'CD2'], 'ILE': 'CD1', 'MET': 'CE',
                      'PHE': 'CZ', 'TRP': 'CH2', 'PRO': 'CG', 'SER': 'OG', 'THR': ['OG1', 'CG2'], 'CYS': 'SG', 'TYR': 'OH',
                      'ASN': ['OD1', 'ND2'], 'GLN': ['OE1', 'NE2'], 'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2'], 'LYS': 'NZ',
                      'ARG': ['NH1', 'NH2'], 'HIS': 'NE2'}

    # https://en.wikipedia.org/wiki/Van_der_Waals_radius
    standard_atom_radius = {'C': 1.7, 'O': 1.52, 'N': 1.55, 'S': 1.8, 'SE': 1.9}
    max_atom_radius = max(standard_atom_radius.values())

    # noinspection PyUnresolvedReferences
    __make_list_of_lists.__func__(standard_residue_atoms)
    # noinspection PyUnresolvedReferences
    __make_list_of_lists.__func__(furthest_atoms)


class AtomPoint():
    def __init__(self, coord, parent, point_id):
        self.coord = coord
        self.parent = parent
        self.id = point_id
        self.deleted = False


class ResidueStatus(Enum):
    modified = 0
    unmodified = 1

    def __str__(self):
        if self.value == 0:
            return 'modified'
        if self.value == 1:
            return 'unmodified'


class SecondaryStructureInfo():
    def __init__(self, chain_id, interval_start, insertion_code_start, interval_end, insertion_code_end,
                 secondary_structure_class):
        self.chain_id = chain_id
        self.interval_start = interval_start
        self.insertion_code_start = insertion_code_start
        self.interval_end = interval_end
        self.insertion_code_end = insertion_code_end
        self.secondary_structure_class = secondary_structure_class

    def __str__(self):
        return 'Interval in {0} from {1}{2} to {3}{4}'.format(self.chain_id, self.interval_start, self.insertion_code_start,
                                                              self.interval_end, self.insertion_code_end)