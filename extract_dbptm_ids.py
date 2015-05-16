import sys
import os.path as path
import collections


# >ACT1_DROME:45:Swiss-Prot 1010711
# PRHQGVmVGMGQK

# extract uniprot ids from dbptm output
def extract_ids():
    filename = sys.argv[1]
    id_set = set()
    with open(filename) as dbptm_output:
        while True:
            line = dbptm_output.readline()
            if not line:
                break
            entry_id = line.split(":", 1)[0][1:]
            id_set.add(entry_id)
            dbptm_output.readline()
    parent_folder = path.split(filename)[0]
    with open(path.join(parent_folder, "result"), "w") as entry_file:
        entry_file.writelines("\n".join(id_set))


# extract uniprot id and stores every modified position in it
def extract_modification_positions():
    filename = sys.argv[1]
    positions = collections.defaultdict(list)
    with open(filename) as dbptm_output:
        while True:
            line = dbptm_output.readline()
            if not line:
                break
            tokens = line.split(":", 2)
            entry_id = tokens[0][1:]
            positions[entry_id].append(tokens[1])
            dbptm_output.readline()
    parent_folder = path.split(filename)[0]
    with open(path.join(parent_folder, "positions"), "w") as entry_file:
        for entry_id, positions in positions.items():
            entry_file.write(entry_id + ":" + ",".join(positions) + "\n")


def main():
    extract_modification_positions()


main()