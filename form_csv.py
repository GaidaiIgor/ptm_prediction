import sys
from os import path
import collections
from Bio import SeqIO


# reads dbptm output with mapping uniprot id -> modified positions, then
# reads fasta file from uniprot with primary structure for each id, then
# knowing modified positions forms csv file with neighbours as predictors


def extract_modification_positions(filename):
    positions = collections.defaultdict(set)
    with open(filename) as dbptm_output:
        while True:
            line = dbptm_output.readline()
            if not line:
                break
            tokens = line.split(":", 2)
            entry_id = tokens[0][1:]
            positions[entry_id].add(int(tokens[1]) - 1)
            dbptm_output.readline()
    return positions


def collect_neighbours(sequence, index, neighbours_amount):
    neighbours = []
    for shift in range(-neighbours_amount, neighbours_amount + 1):
        if shift == 0:
            continue
        neighbour_index = index + shift
        if neighbour_index < 0 or neighbour_index >= len(sequence):
            neighbours.append("-")
            continue
        neighbours.append(sequence[index + shift])
    return neighbours


def main():
    positions_filename = sys.argv[1]
    sequences_filename = sys.argv[2]
    neighbours_amount = 3
    modification_positions = extract_modification_positions(positions_filename)
    all_features = []
    with open(sequences_filename) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            entry_id = record.name.split("|")[-1]
            for i in range(len(record.seq)):
                amino_acid = record.seq[i]
                if amino_acid != "M":
                    continue
                features = [entry_id]
                if i in modification_positions[entry_id]:
                    features.append("modified")
                else:
                    features.append("unmodified")
                for neighbour in collect_neighbours(record.seq, i, neighbours_amount):
                    features.append(neighbour)
                all_features.append(features)
    # parent_folder = path.split(sequences_filename)[0]
    header = ["filename", "status"] + list(map(lambda i: "type" + str(i), range(neighbours_amount * 2)))
    with open("/home/dartlenin/R_Projects/ptm_prediction/data_primary_dbptm.csv", "w+") as output:
        output.write(",".join(header) + "\n")
        for observation in all_features:
            output.write(",".join(observation) + "\n")


main()