import sys
import time

# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import Axes3D
from FeatureCollector import FeatureCollector


def square_distance(point1, point2):
    return sum([(coord1 - coord2) ** 2 for (coord1, coord2) in zip(point1, point2)])


def test_residue_sas_points(residue, structure, solvent_radius):
    for atom in residue.atoms:
        for point in atom.sas_points:
            for other_atom in structure.get_atoms():
                if other_atom != atom:
                    if square_distance(point.coord, other_atom.coord) < (other_atom.radius + solvent_radius) ** 2:
                        return False
    return True


def main():
    data_path = sys.argv[1]
    with open('output1.csv', 'w') as output:
        feature_collector = FeatureCollector()
        feature_collector.export_features(output, data_path)


# check SNP: 3NBJ - 634 A
start = time.clock()
# a = [1, [2, 3], [1, [2, 3, [1, [2, 3]]]]]
# print(list(FeatureCollector.flat_list(a)))
# b(s = 3, d = 1)
# c()
main()
# cProfile.run('main()', '1')
print('{0}s elapsed'.format(time.clock() - start))