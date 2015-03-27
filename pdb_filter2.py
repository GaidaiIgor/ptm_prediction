import sys
import os
import os.path as path
import shutil


def modification_names(pdb_file, residue_name):
    modres_names = set()
    directories = set()
    for line in pdb_file:
        record_type = line[:6].strip()
        # collect all modres short names
        if record_type == 'MODRES':
            std_res_name = line[24:27].strip()
            if std_res_name == residue_name:
                modres_name = line[12:15].strip()
                modres_names.add(modres_name)
        if record_type == 'HETNAM':
            het_name = line[11:14].strip()
            # if it's a modified residue name
            if het_name in modres_names:
                het_full_name = line[15:].strip()
                directories.add(het_full_name)
    return directories


def main():
    print_step = 100
    data_path = sys.argv[1]
    residue_name = sys.argv[2]
    filenames = os.listdir(data_path)

    # filenames = filenames[:1]
    # filenames = ['132L.pdb']

    files_proceeded = 0
    for filename in filenames:
        next_file_path = path.join(data_path, filename)
        with open(next_file_path) as next_file:
            directories = modification_names(next_file, residue_name)
        for directory in directories:
            full_destination_path = path.join(data_path, directory)
            if not path.exists(full_destination_path):
                os.mkdir(full_destination_path)
            shutil.copy(next_file_path, full_destination_path)
        files_proceeded += 1
        if files_proceeded % print_step == 0:
            print('%d files out of %d files proceeded' % (files_proceeded, len(filenames)))


main()