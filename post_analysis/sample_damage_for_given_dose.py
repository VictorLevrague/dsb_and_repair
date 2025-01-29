import random
import re
from collections import defaultdict
import pandas as pd
from scipy.stats import poisson
import numpy as np
import time
import os.path

# PATH = "/sps/gdrmi2b/levrague/dsb_and_repair"
PATH = "/home/levrague/Documents/PostDoc/dsb_and_repair"
IRRADIATION_ID = 1 #1 (HC nucleus) or 2 (decompacted nucleus)
NB_TRACKS_1_GY_DICT = {"Cs_137": 741592, "H_250": 4390, "H_150": 3166, "He_250": 1110, "C_290": 135,"O_350": 84, "O_55": 24, "Si_170": 17, "Ti_300": 10, "Fe_600": 11,
                       "Fe_450": 9, "Fe_300": 7}
NB_SIMULATED_TRACKS = {"Cs_137": 738000, "H_250": 4142, "H_150": 2976,"He_250": 1035, "C_290": 125,"O_350": 78,"O_55": 22, "Si_170": 17, "Ti_300": 10,
                       "Fe_600": 9, "Fe_450": 8, "Fe_300": 7}
ION_LET = {"Cs_137":0.2, "H_250": 0.39, "H_150": 0.54, "He_250": 1.54, "C_290": 12.74,
           "O_350": 20.47, "O_55": 72.84, "Si_170": 96.6, "Ti_300": 168.06,
           "Fe_600": 170.62, "Fe_450": 191.57, "Fe_300": 234.71} #keV/µm
SIMULATED_DOSE = 1 #Gy
SAMPLE_DOSES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1] #Gy
NB_SIMULATIONS = 1000
NB_SIMULATION_SAMPLES = 10000
SEED = 40
# "Cs_137_EC", "Cs_137_HC"
IONS = ["H_250MeV_EC", "H_250MeV_HC",
        "He_250MeV_EC", "He_250MeV_HC", "O_350MeV_EC", "O_350MeV_HC",
        "O_55MeV_EC", "O_55MeV_HC", "Si_170MeV_EC", "Si_170MeV_HC",
        "Ti_300MeV_EC", "Ti_300MeV_HC", "Fe_600MeV_EC", "Fe_600MeV_HC",
        "Fe_450MeV_EC", "Fe_450MeV_HC", "Fe_300MeV_EC", "Fe_300MeV_HC"]

def nb_ions_poisson(dose, ion):
    # AREA = 16 * 16 #µm²
    # UNIT_FACTOR = 6.24 #keV.kg.J^-1.µm^-3
    # lambda_poisson = UNIT_FACTOR * dose * AREA / ION_LET[ion]
    lambda_poisson = int(NB_TRACKS_1_GY_DICT[ion] * dose)
    return poisson.rvs(lambda_poisson)

def compute_nb_damage(damage_file_dict, ion, dose, all_possible_tracks):
    sample_damage_dict = defaultdict(list)
    damage_dict = defaultdict()
    damage_dict_error = defaultdict()
    for nb_sample in range(NB_SIMULATION_SAMPLES):
        nb_tracks_to_sample = nb_ions_poisson(dose, normalize_key(ion))
        sampled_tracks = random.sample(all_possible_tracks, nb_tracks_to_sample)
        for damage_modifier in damage_modifiers:
            sample_damage_dict[damage_modifier].append(
                np.sum(damage_file_dict[track][damage_modifier] for track in sampled_tracks))
    for damage_modifier in damage_modifiers:
        sample_damages = np.array(sample_damage_dict[damage_modifier])
        damage_dict[damage_modifier] = np.mean(sample_damages)
        damage_dict_error[damage_modifier] = np.std(sample_damages) / np.sqrt(NB_SIMULATION_SAMPLES)
    return damage_dict, damage_dict_error

def convert_damage_file_in_dict(damage_type, ion):
    damage_file_dict = defaultdict(lambda: defaultdict(int))
    damage_file_name = f"{PATH}/post_analysis/analysis_output/{ion}/Irradiation{IRRADIATION_ID}/{SIMULATED_DOSE}Gy/List_{damage_type}.txt"
    with open(damage_file_name, "r") as damage_file:
        [next(damage_file) for _ in range(0, 3)]
        for line in damage_file:
            if not line.strip():  #remove white spaces
                continue
            parts = line.split()
            run_id = parts[0]
            event_id = parts[1]
            damage_file_dict[(int(run_id), int(event_id))]["Total"] += 1
            if int(parts[8]) > 1:
                damage_file_dict[(int(run_id), int(event_id))]["Complex"] += 1
            if int(parts[9]):
                damage_file_dict[(int(run_id), int(event_id))]["Direct"] += 1
            if int(parts[10]):
                damage_file_dict[(int(run_id), int(event_id))]["Indirect"] += 1
    return damage_file_dict

def normalize_key(key):
    match = re.match(r"^[A-Za-z]+_\d+", key)
    return match.group(0) if match else None

def set_rng():
    np.random.seed(SEED)
    random.seed(SEED)

if __name__ == '__main__':
    set_rng()
    output_dict = defaultdict(list)
    start = time.time()
    damage_modifiers = ["Total", "Direct", "Indirect", "Complex"]
    for ion in IONS:
        all_possible_tracks = [(sim_id, track_id) for sim_id in range(NB_SIMULATIONS)
                               for track_id in range(NB_SIMULATED_TRACKS[normalize_key(ion)])]
        damage_file_dict_ssb = convert_damage_file_in_dict("SSB", ion)
        damage_file_dict_dsb = convert_damage_file_in_dict("DSB", ion)
        for dose in SAMPLE_DOSES:
            ssb_dict, ssb_dict_error = compute_nb_damage(damage_file_dict_ssb, ion, dose, all_possible_tracks)
            dsb_dict, dsb_dict_error = compute_nb_damage(damage_file_dict_dsb, ion, dose, all_possible_tracks)
            output_dict["Ion"].append(ion)
            output_dict["LET"].append(ION_LET[normalize_key(ion)])
            output_dict["Dose"].append(dose)
            for damage_modifier in damage_modifiers:
                output_dict[f"SSB_{damage_modifier}"].append(ssb_dict[damage_modifier])
                output_dict[f"SSB_{damage_modifier}_error"].append(ssb_dict_error[damage_modifier])
                output_dict[f"DSB_{damage_modifier}"].append(dsb_dict[damage_modifier])
                output_dict[f"DSB_{damage_modifier}_error"].append(dsb_dict_error[damage_modifier])
            print("dose: ", dose, "analyzed")
        print(f"{ion} is analysed\n")
    output_df = pd.DataFrame(output_dict)
    output_df.to_csv("damage_sample_analysis.csv", index=False)
    print("end: ", time.time() - start)