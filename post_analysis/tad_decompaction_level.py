from collections import defaultdict
import numpy as np
import os

nb_simulations = 1
ion = "O_350_MeV"
total_nb_base_pairs = 6101655710

def convert_tad_base_pair_file_to_dictionary():
    tad_bp_conversion_dict = defaultdict(int)
    with open("tad_base_pairs_hc.txt", "r") as tad_bp_file:
        for line in tad_bp_file:
            parts = line.split()
            if (len(parts) >= 2):
                tad_bp_conversion_dict[int(parts[0])] = int(parts[1])
    return tad_bp_conversion_dict

def compute_nb_base_pair_decompacted(tad_hit_filename):
    tad_bp_conversion_dict = convert_tad_base_pair_file_to_dictionary()
    np_base_pairs_decompacted = 0
    with open(tad_hit_filename, "r") as tad_hit_file:
        for line in tad_hit_file:
            tad_id = line.split()[0]
            np_base_pairs_decompacted += tad_bp_conversion_dict[int(float(tad_id))]
    return np_base_pairs_decompacted

def compute_decompaction_fraction(tad_hit_filename):
    return compute_nb_base_pair_decompacted(tad_hit_filename) / total_nb_base_pairs

def compute_average_decompaction_fraction(damage_type):
    decompaction_fraction_array = np.array([])
    for copy_nb in range (0, nb_simulations):
        # print(f"Reading tad hit file n°{copy_nb}")
        file_name = f"tad_hits/{ion}/tad_hits_{copy_nb}_{damage_type}.txt"
        decompaction_fraction_array = np.append(decompaction_fraction_array, compute_decompaction_fraction(file_name))
    return np.mean(decompaction_fraction_array)

def compute_nb_tad_hits(tad_hit_filename):
    with open(tad_hit_filename, "r") as tad_hit_file:
        lines = tad_hit_file.readlines()
    return len(lines)

def compute_mean_nb_tad_hits(damage_type):
    tad_hits_array = np.array([])
    for copy_nb in range (0, nb_simulations):
        file_name = f"tad_hits/{ion}/tad_hits_{copy_nb}_{damage_type}.txt"
        tad_hits_array = np.append(tad_hits_array, compute_nb_tad_hits(file_name))
    return np.mean(tad_hits_array)

if __name__ == '__main__':
    damage_types = ["SSB", "DSB", "ALL"]
    for damage_type in damage_types:
        mean_decompaction_fraction = compute_average_decompaction_fraction(damage_type)
        print(f"Mean genome decompaction percentage for {ion} with {damage_type} is:  {round(mean_decompaction_fraction * 100, 2)} %")
        print(f"Mean nb of decompacted DNA for {ion} with {damage_type} is:  {round(mean_decompaction_fraction * 6101655710 / 1e6, 2) } Mbp")
        print(f"Mean nb of tad hits for {ion} with {damage_type} is:  {round(compute_mean_nb_tad_hits(damage_type), 2)} TADs")
        print()
