import sys
import os
import os.path as path

import Bio.PDB as PDB


def get_unmodified_residue_codes(pdb_file):
    modres_names = {}
    for line in pdb_file:
        if line[:6].strip() == 'MODRES':
            modres_name = line[12:15].strip()
            std_residue_name = line[24:27].strip()
            modres_names[modres_name] = std_residue_name
    return modres_names


def process_directory(data_path):
    print_step = 10

    pdb_parser = PDB.PDBParser(PERMISSIVE = False, QUIET = True)
    filenames = os.listdir(data_path)

    # filenames = filenames[:1]
    # filenames = ['1E93.pdb']

    files_proceeded = 0
    for filename in filenames:
        next_file_path = path.join(data_path, filename)
        # if path.isdir(next_file_path):
        # files_proceeded += process_directory(next_file_path)
        if path.islink(next_file_path):
            process_directory(os.readlink(next_file_path))
        else:
            with open(next_file_path) as next_file:
                unmod_res_codes = get_unmodified_residue_codes(next_file)
                next_file.seek(0)
                structure = pdb_parser.get_structure(filename, next_file)
            for residue in structure.get_residues():
                # we don't need water
                if residue.id[0] == 'W':
                    continue
                # if it's a modified residue
                # so we should build spheres only around standard atoms
                # because we want to predict modification that didn't happen yet
                elif residue.id[0].startswith('H_'):
                    # skip non amino acid het atoms
                    if residue.resname not in unmod_res_codes:
                        continue
                    for possible_atom_names in standard_residue_atoms[unmod_res_codes[residue.resname]]:
                        # create 1-element tuple to avoid iterating over string letters
                        if not isinstance(possible_atom_names, tuple):
                            possible_atom_names = (possible_atom_names, )
                        for atom_name in possible_atom_names:
                            if residue.has_id(atom_name):
                                break
                        # residue lacks some necessary atom
                        else:
                            print('file: ' + structure.id + ', chain: ' + residue.parent.id +
                                  ', residue: ' + str(residue.id) + ' lacks some necessary atom: ' + str(possible_atom_names))
                # if residue is not modified then just iterate over all atoms and check that every atom present in standard atoms
                else:
                    for atom in residue:
                        # skip hydrogen atoms
                        if atom.element == 'H':
                            continue
                        if atom.id == 'OXT':
                            continue
                        # check if atom is in standard residue atoms
                        for possible_atom_names in standard_residue_atoms[residue.resname]:
                            if not isinstance(possible_atom_names, tuple):
                                possible_atom_names = (possible_atom_names, )
                            if atom.id in possible_atom_names:
                                break
                        else:
                            print('file: ' + structure.id + ', chain: ' + residue.parent.id +
                                  ', residue: ' + str(residue.id) + ', atom: ' + atom.id + ' is not in standard atoms')

            files_proceeded += 1
        if files_proceeded % print_step == 0:
            print('{0} out of {1} files in {2} proceeded'.format(files_proceeded, len(filenames), data_path))
    return files_proceeded


def main():
    data_path = sys.argv[1]
    process_directory(data_path)


standard_backbone_atoms = ['N', 'CA', 'C', ('O', 'S')]
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
# check SNP: 3NBJ - 634 A
# check multiple models
main()