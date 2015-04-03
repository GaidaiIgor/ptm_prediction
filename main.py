import sys
import time

import matplotlib.pyplot as pyplot


# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import Axes3D
from Structure import ProteinStructure
from FeatureCollector import FeatureCollector


def get_surface_around(point, surface_kdtree, ):
    surface_around = []


def get_near_surface_stats(residue, surface_kdtree, query_radius, tree_depth):
    furthest_atom = ProteinStructure.get_furthest_atom(residue)
    if furthest_atom is None:
        return ['']
    nearest_point_info = surface_kdtree.query(furthest_atom.coord)
    is_used = [False] * len(surface_kdtree.data)


def plot_points(points):
    figure = pyplot.figure()
    axes = figure.add_subplot(111, projection = '3d')
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    zs = [point[2] for point in points]
    axes.scatter(xs, ys, zs)
    pyplot.show()


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


def a(s = 2, d = 3):
    print(s, d)

def b(s = 1, **kwargs):
    print(s)
    a(s, **kwargs)


# check SNP: 3NBJ - 634 A
# check multiple models
start = time.clock()
# b(s = 3, d = 1)
main()
# cProfile.run('main()', '5')
print('{0} time elapsed'.format(time.clock() - start))