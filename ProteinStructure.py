from enum import Enum
import itertools
import os

from scipy.spatial import KDTree
from Bio import PDB


class ProteinStructure(object):
    def __init__(self, pdb_file=None, solvent_radius=1.5, surface_triangulation_density=1):
        # self.atom_points = []
        self.atoms = []
        self.edge_coords = []
        self.filename = ''
        self.mod_code_to_unmod_code = {}
        self.mod_name_to_mod_code = {}
        self.protein_residues = []
        self.structure = None
        self.ses = []
        self.surface_kdtree = None
        self.atom_kdtree = None
        self.residues_key_to_pos = {}
        self.surface_triangulation_density = surface_triangulation_density
        self.solvent_radius = solvent_radius
        if pdb_file is not None:
            self.parse_structure(pdb_file)

    def get_surface(self):
        return self.ses

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

    def get_residue_by_key(self, key):
        return self.get_residues()[self.residues_key_to_pos[key]]

    def parse_msms_vertices_file(self, vertices_file):
        atoms = self.get_atoms()
        for line in vertices_file:
            line_content = line.split()
            point_coord = [float(token) for token in line_content[0:3]]
            point_normal_vector = [float(token) for token in line_content[3:6]]
            nearest_atom_index = int(line_content[7]) - 1
            point_nearest_atom = atoms[nearest_atom_index]
            surface_point = SurfacePoint(point_coord, point_normal_vector, point_nearest_atom)
            self.ses.append(surface_point)

    def parse_atom_area(self, atom_area_file):
        atoms = self.get_atoms()
        i = 0
        _ = next(atom_area_file)
        for line in atom_area_file:
            line_content = [float(token) for token in line.split()]
            next_atom = atoms[i]
            next_atom.ses = line_content[1]
            next_atom.sas = line_content[2]
            i += 1
        assert i == len(atoms), 'Area is defined not for each atom'

    def run_msms(self):
        with open('xyzr', 'w') as xyzr_file:
            self.to_xyzr(xyzr_file)
        surface_path_prefix = 'surface'
        atom_area_path_prefix = 'atom_area'
        msms_command = 'msms -no_header -probe_radius {} -density {} -noh -if {} -of {} -af {} > {}'. \
            format(self.solvent_radius, self.surface_triangulation_density, xyzr_file.name, surface_path_prefix,
                   atom_area_path_prefix, os.devnull)
        os.system(msms_command)
        msms_produced_files = [surface_path_prefix + '.vert', surface_path_prefix + '.face', atom_area_path_prefix + '.area']
        surface_vertices_file_path = msms_produced_files[0]
        atom_area_file_path = msms_produced_files[2]

        with open(surface_vertices_file_path) as surface_vertices_file:
            self.parse_msms_vertices_file(surface_vertices_file)
        with open(atom_area_file_path) as atom_area_file:
            self.parse_atom_area(atom_area_file)

        os.remove(xyzr_file.name)
        for file_path in msms_produced_files:
            os.remove(file_path)

    def build_ses(self):
        self.run_msms()

    @staticmethod
    def atom_to_xyzr(atom):
        return '{} {} {} {}'.format(atom.coord[0], atom.coord[1], atom.coord[2], atom.radius)

    def to_xyzr(self, xyzr_file):
        atoms = self.get_atoms()
        for atom in atoms:
            xyzr_file.write(ProteinStructure.atom_to_xyzr(atom) + '\n')
        xyzr_file.flush()

    def __add_atom_radii(self):
        for atom in self.get_atoms():
            if atom.element in ProteinStructure.standard_atom_radius:
                atom.radius = ProteinStructure.standard_atom_radius[atom.element]
            # shouldn't happen
            else:
                print('{0} at {1} in {2}: atom radius is unknown'.format(atom, atom.parent.id[1], self.filename))

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

        all_residues = self.get_residues()
        # now let's change it for specified intervals
        for info in secondary_structure_info:
            start_key = ProteinStructure.make_residue_key(info.chain_id, info.interval_start, info.insertion_code_start)
            end_key = ProteinStructure.make_residue_key(info.chain_id, info.interval_end, info.insertion_code_end)
            pos_start = self.residues_key_to_pos[start_key]
            pos_end = self.residues_key_to_pos[end_key]
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

    def get_surface_points(self):
        return self.ses

    def get_surface_kdtree(self):
        if self.surface_kdtree is None:
            surface_points = self.get_surface_points()
            points = [point.coord for point in surface_points]
            self.surface_kdtree = KDTree(points)
        return self.surface_kdtree

    def __init_key_to_pos(self):
        all_residues = self.get_residues()
        for i in range(len(all_residues)):
            self.residues_key_to_pos[ProteinStructure.__residue_to_key(all_residues[i])] = i

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
            self.protein_residues.append(residue)
        self.__init_key_to_pos()

    @staticmethod
    def __make_list_of_lists(target_dict: dict):
        for key, value in target_dict.items():
            if not isinstance(value, list):
                target_dict[key] = [value]

    @staticmethod
    def make_residue_key(residue_chain, residue_number, insertion_code=' '):
        return residue_chain, residue_number, insertion_code

    @staticmethod
    def __residue_to_key(residue):
        return residue.parent.id, residue.id[1], residue.id[2]

    beta_sheet_class_shift = 12
    no_secondary_structure_class = 14
    # solvent radius = 1.2 + 0.96 where 1.2 is oh bond length in water molecule and 0.96 is van der waals radius for hydrogen
    # https://ru.wikipedia.org/wiki/%D0%92%D0%BE%D0%B4%D0%B0#/media/File:Water_molecule_dimensions.svg
    # https://en.wikipedia.org/wiki/Van_der_Waals_radius
    # default_solvent_radius = 2.16

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


class SurfacePoint():
    def __init__(self, coord, normal_vector, closest_atom):
        self.coord = coord
        self.normal_vector = normal_vector
        self.closest_atom = closest_atom


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