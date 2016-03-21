import sys
import time

# noinspection PyUnresolvedReferences
# from mpl_toolkits.mplot3d import Axes3D
from FeatureCollector import FeatureCollector


def main():
    neighbours_amount = 6
    solvent_radius = 1.5
    data_path = sys.argv[1]
    output_file_path = r"E:\Dropbox\masters_work\r_project\met6_1.5.csv"
    modification_full_name = sys.argv[2]
    with open(output_file_path, "w") as output:
        feature_collector = FeatureCollector(neighbours_amount=neighbours_amount,
                                             protein_structure_args={"solvent_radius": solvent_radius})
        feature_collector.export_features(output, data_path, modification_full_name)
    # with open(path.join(output_path, "met_primary6.csv"), "w") as output:
    #     feature_collector = FeatureCollector(neighbours_amount=neighbours_amount, use_tertiary=False)
    #     feature_collector.export_features(output, data_path)


start = time.clock()
main()
print('{0}s elapsed'.format(time.clock() - start))
