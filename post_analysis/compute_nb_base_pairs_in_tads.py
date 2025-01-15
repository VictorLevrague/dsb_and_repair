from collections import defaultdict
import numpy as np
import os

filename_hc = "dnafabric_geometries/geometrie_fibroblast_hc_floriane.fab2g4dna"
filename_ec = "dnafabric_geometries/geometrie_fibroblast_ec_floriane.fab2g4dna"

bp_per_voxel_dict = {"VoxelStraight": 3594,
                     "VoxelLeft": 2402,
                     "VoxelRight": 2404,
                     "VoxelUp": 2408,
                     "VoxelDown": 2415,
                     "VoxelStraight2": 2011,
                     "VoxelLeft2": 1646,
                     "VoxelRight2": 1646,
                     "VoxelUp2": 1637,
                     "VoxelDown2": 1660
                     }

def write_tad_bp_file(filename):
    with open(filename, "r") as geometry_file:
        _ = [geometry_file.readline() for _ in range(11)]
        tad_bp_dict = compute_nb_bd_in_tad(geometry_file)
    return tad_bp_dict

def compute_total_nb_base_pairs(tad_bp_dict):
    nb_base_pairs = 0
    for key in tad_bp_dict.keys():
        nb_base_pairs += tad_bp_dict[key]
    return nb_base_pairs

def compute_nb_bd_in_tad(geometry_file):
    tad_bp_dict = defaultdict(int)
    for line in geometry_file:
        parts = line.split()
        if (len(parts) >= 4):
            voxel_type = parts[1]
            tad_id = parts[3]
            tad_bp_dict[tad_id] += bp_per_voxel_dict[voxel_type]
    return tad_bp_dict

if __name__ == '__main__':
    bp_hc_dict = write_tad_bp_file(filename_hc)
    print("Number of base pairs in hc geometry: ", compute_total_nb_base_pairs(bp_hc_dict))
    bp_ec_dict = write_tad_bp_file(filename_ec)
    print("Number of base pairs in ec geometry: ", compute_total_nb_base_pairs(bp_ec_dict))
    #
    with open("tad_base_pairs_hc.txt", 'w') as file:
        for tad_id in bp_hc_dict.keys():
            file.write(f"{tad_id}\t{bp_hc_dict[tad_id]}\n")