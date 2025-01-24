from collections import defaultdict
import pandas as pd
import numpy as np
import re
import os

# PATH = "/sps/gdrmi2b/levrague/dsb_and_repair"
PATH = "/home/levrague/Documents/PostDoc/dsb_and_repair"
NB_SIMULATIONS = 1000
IONS = ["H_250", "He_250", "O_55", "O_350", "Si_170", "Ti_300",
        "Fe_300", "Fe_450", "Fe_600"]
ION_LET = {"H_250": 0.39, "H_150": 0.54, "He_250": 1.54, "C_290": 12.74,
           "O_350": 20.47, "O_55": 72.84, "Si_170": 96.6, "Ti_300": 168.06,
           "Fe_600": 170.62, "Fe_450": 191.57, "Fe_300": 234.71} #keV/µm
DOSES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1] #Gy
# TOTAL_NB_BASE_PAIRS = 6101655710

def convert_tad_base_pair_file_to_dictionary():
    tad_bp_conversion_dict = defaultdict(int)
    with open("tad_base_pairs_hc.txt", "r") as tad_bp_file:
        for line in tad_bp_file:
            parts = line.split()
            if (len(parts) >= 2):
                tad_bp_conversion_dict[int(parts[0])] = int(parts[1])
    return tad_bp_conversion_dict

def compute_nb_base_pair_decompacted_and_tad_hit(tad_hit_filename):
    tad_bp_conversion_dict = convert_tad_base_pair_file_to_dictionary()
    np_base_pairs_decompacted = 0
    nb_tad_hit = 0
    with open(tad_hit_filename, "r") as tad_hit_file:
        for line in tad_hit_file:
            tad_id = line.split()[0]
            np_base_pairs_decompacted += tad_bp_conversion_dict[int(float(tad_id))]
            nb_tad_hit += 1
    return np_base_pairs_decompacted, nb_tad_hit

# def compute_decompaction_fraction(tad_hit_filename):
#     return compute_nb_base_pair_decompacted(tad_hit_filename) / TOTAL_NB_BASE_PAIRS
#
# def compute_average_decompaction_fraction(damage_type):
#     decompaction_fraction_array = np.array([])
#     for copy_nb in range (0, NB_SIMULATIONS):
#         print(f"Reading tad hit file n°{copy_nb}")
#         file_name = f"{output_folder}/tad_hits_{copy_nb}_{damage_type}.txt"
#         decompaction_fraction_array = np.append(decompaction_fraction_array, compute_decompaction_fraction(file_name))
#     return np.mean(decompaction_fraction_array)

def generate_array_decompacted_dna_and_nb_tad_hits(ion, dose, damage_type):
    output_folder = f"{PATH}/post_analysis/tad_hits/{ion}/{dose}Gy"
    decompacted_dna_array = np.array([])
    nb_tad_hit_array = np.array([])
    for copy_nb in range(0, NB_SIMULATIONS):
        print(f"Reading tad hit file n°{copy_nb}")
        file_name = f"{output_folder}/tad_hits_{copy_nb}_{damage_type}.txt"
        nb_bp_decompacted, nb_tad_hit = compute_nb_base_pair_decompacted_and_tad_hit(file_name)
        decompacted_dna_array = np.append(decompacted_dna_array, nb_bp_decompacted)
        nb_tad_hit_array = np.append(nb_tad_hit_array, nb_tad_hit)
    return decompacted_dna_array, nb_tad_hit_array

def compute_nb_tad_hits(tad_hit_filename):
    with open(tad_hit_filename, "r") as tad_hit_file:
        lines = tad_hit_file.readlines()
    return len(lines)

def normalize_key(key):
    match = re.match(r"^[A-Za-z]+_\d+", key)
    return match.group(0) if match else None

# def compute_mean_nb_tad_hits(damage_type):
#     tad_hits_array = np.array([])
#     for copy_nb in range (0, NB_SIMULATIONS):
#         file_name = f"{output_folder}/tad_hits_{copy_nb}_{damage_type}.txt"
#         tad_hits_array = np.append(tad_hits_array, compute_nb_tad_hits(file_name))
#     return np.mean(tad_hits_array)

if __name__ == '__main__':
    output_dict = defaultdict(lambda: np.array([]))
    damage_types = ["SSB", "DSB", "ALL"]
    for ion in IONS:
        for dose in DOSES:
            decompacted_dna_ssb, tad_hit_ssb = generate_array_decompacted_dna_and_nb_tad_hits(ion, dose, "SSB")
            decompacted_dna_dsb, tad_hit_dsb = generate_array_decompacted_dna_and_nb_tad_hits(ion, dose, "DSB")
            decompacted_dna_all, tad_hit_all = generate_array_decompacted_dna_and_nb_tad_hits(ion, dose, "ALL")
            output_dict["Ion"] = np.append(output_dict["Ion"], ion)
            output_dict["LET"] = np.append(output_dict["LET"], ION_LET[normalize_key(ion)])
            output_dict["Dose"] = np.append(output_dict["Dose"], dose)
            output_dict["DecompDNA_SSB"] = np.append(output_dict["DecompDNA_SSB"], np.mean(decompacted_dna_ssb))
            output_dict["DecompDNA_SSB_error"] = np.append(output_dict["DecompDNA_SSB_error"], np.std(decompacted_dna_ssb)/np.sqrt(NB_SIMULATIONS))
            output_dict["DecompDNA_DSB"] = np.append(output_dict["DecompDNA_DSB"], np.mean(decompacted_dna_dsb))
            output_dict["DecompDNA_DSB_error"] = np.append(output_dict["DecompDNA_DSB_error"], np.std(decompacted_dna_dsb)/np.sqrt(NB_SIMULATIONS))
            output_dict["DecompDNA_ALL"] = np.append(output_dict["DecompDNA_ALL"], np.mean(decompacted_dna_all))
            output_dict["DecompDNA_ALL_error"] = np.append(output_dict["DecompDNA_ALL_error"],
                                                           np.std(decompacted_dna_all) / np.sqrt(NB_SIMULATIONS))
            output_dict["NbTadHits_SSB"] = np.append(output_dict["NbTadHits_SSB"], np.mean(tad_hit_ssb))
            output_dict["NbTadHits_SSB_error"] = np.append(output_dict["NbTadHits_SSB_error"], np.std(tad_hit_ssb)/np.sqrt(NB_SIMULATIONS))
            output_dict["NbTadHits_DSB"] = np.append(output_dict["NbTadHits_DSB"], np.mean(tad_hit_dsb))
            output_dict["NbTadHits_DSB_error"] = np.append(output_dict["NbTadHits_DSB_error"], np.std(tad_hit_dsb)/np.sqrt(NB_SIMULATIONS))
            output_dict["NbTadHits_ALL"] = np.append(output_dict["NbTadHits_ALL"], np.mean(tad_hit_all))
            output_dict["NbTadHits_ALL_error"] = np.append(output_dict["NbTadHits_ALL_error"],
                                                           np.std(tad_hit_all) / np.sqrt(NB_SIMULATIONS))
            print(f"{dose} is analysed")
        print(f"{ion} is analysed\n")
    output_df = pd.DataFrame(output_dict)
    output_df.to_csv("tad_hit_analysis.csv", index=False)