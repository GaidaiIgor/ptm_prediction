import sys
from os import path
# reads uniprot mapping file like
# From	To
# PTPRG_MOUSE	3JXG
# PTPRG_MOUSE	3KLD
# NAP1_YEAST	2AYU
# NAP1_YEAST	2Z2R
# and then select first pdb_id for each uniprot id and prints it


def main():
    filename = sys.argv[1]
    result = {}
    with open(filename) as pdb_id_map:
        for line in pdb_id_map:
            line_tokens = line.split("\t")
            result[line_tokens[0]] = line_tokens[1]
    parent = path.split(filename)[0]
    with open(path.join(parent, "pdb_id"), "w") as output:
        for pdb_id in result.values():
            output.write(pdb_id)

main()
