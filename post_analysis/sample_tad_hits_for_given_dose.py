import random
import re
from collections import defaultdict
from scipy.stats import poisson
import numpy as np
import os.path

# PATH = "/sps/gdrmi2b/levrague/dsb_and_repair"
PATH = "/home/levrague/Documents/PostDoc/dsb_and_repair"
IRRADIATION_ID = 1 #1 (HC nucleus) or 2 (decompacted nucleus)
NB_TRACKS_1_GY_DICT = {"Cs-137": 741592,"H_250": 4390, "H_150": 3166, "He_250": 1110, "C_290": 135,"O_350": 84, "O_55": 24, "Si_170": 17, "Ti_300": 10, "Fe_600": 11,
                       "Fe_450": 9, "Fe_300": 7}
NB_SIMULATED_TRACKS = {"Cs-137": 738000, "H_250": 4142, "H_150": 2976,"He_250": 1035, "C_290": 125,"O_350": 78,"O_55": 22, "Si_170": 17, "Ti_300": 10,
                       "Fe_600": 9, "Fe_450": 8, "Fe_300": 7}
ION_LET = {"Cs-137":0.2,"H_250": 0.39, "H_150": 0.54, "He_250": 1.54, "C_290": 12.74,
           "O_350": 20.47, "O_55": 72.84, "Si_170": 96.6, "Ti_300": 168.06,
           "Fe_600": 170.62, "Fe_450": 191.57, "Fe_300": 234.71} #keV/µm
SIMULATED_DOSE = 1 #Gy
# SAMPLE_DOSES = [1] #Gy
SAMPLE_DOSES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1] #Gy
NB_SIMULATIONS = 1000
NB_SIMULATION_SAMPLES = 1000
SEED = 42
# "Cs-137_HC"
IONS = ["H_250MeV_HC", "He_250MeV_HC", "O_350MeV_HC", "O_55MeV_HC","Si_170MeV_HC",
        "Ti_300MeV_HC", "Fe_600MeV_HC", "Fe_450MeV_HC", "Fe_300MeV_HC"]

def nb_ions_poisson(dose, ion):
    # AREA = 16 * 16 #µm²
    # UNIT_FACTOR = 6.24 #keV.kg.J^-1.µm^-3
    # lambda_poisson = UNIT_FACTOR * dose * AREA / ION_LET[ion]
    lambda_poisson = int(NB_TRACKS_1_GY_DICT[ion] * dose)
    return poisson.rvs(lambda_poisson)

def write_tad_hit_files_with_sampling(ion):
    damage_file_dict_ssb = convert_damage_file_in_tad_id_dict("SSB", ion)
    damage_file_dict_dsb = convert_damage_file_in_tad_id_dict("DSB", ion)
    for dose in SAMPLE_DOSES:
        output_folder = f"{PATH}/post_analysis/tad_hits/{normalize_key(ion)}/{dose}Gy"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for nb_sample in range(0, NB_SIMULATION_SAMPLES):
            nb_tracks_to_sample = nb_ions_poisson(dose, normalize_key(ion))
            tad_hits_sample = compute_tad_hit_array_for_one_sample(damage_file_dict_ssb, damage_file_dict_dsb, nb_tracks_to_sample)
            write_tad_hit_file("SSB", np.unique(tad_hits_sample["SSB"]), output_folder, nb_sample)
            write_tad_hit_file("DSB", np.unique(tad_hits_sample["DSB"]), output_folder, nb_sample)
            write_tad_hit_file("ALL", np.unique(np.concatenate((tad_hits_sample["SSB"], tad_hits_sample["DSB"]))), output_folder, nb_sample)
            if (nb_sample + 1) % 100 == 0:
                print(f"Sample n°{nb_sample + 1} / {NB_SIMULATION_SAMPLES}")
        print(f"{dose} is analysed")
    return

def compute_tad_hit_array_for_one_sample(damage_file_dict_ssb, damage_file_dict_dsb, nb_tracks_to_sample):
    tad_hits_sample_dict = {"SSB": np.array([]), "DSB": np.array([]), "ALL": np.array([])}
    already_sampled_tracks = defaultdict(int)  # Simulation_id: track_id
    for nb_track in range(0, nb_tracks_to_sample):
        simulation_id = random.randint(0, NB_SIMULATIONS - 1)
        track_id = random.randint(0, NB_TRACKS_1_GY_DICT[normalize_key(ion)])
        while (simulation_id, track_id) in already_sampled_tracks.items():
            simulation_id = random.randint(0, NB_SIMULATIONS - 1)
            track_id = random.randint(0, NB_TRACKS_1_GY_DICT[normalize_key(ion)])
        already_sampled_tracks[simulation_id] = track_id
        tad_hits_sample_dict["SSB"] = np.append(tad_hits_sample_dict["SSB"], damage_file_dict_ssb[(simulation_id, track_id)])
        tad_hits_sample_dict["DSB"] = np.append(tad_hits_sample_dict["DSB"], damage_file_dict_dsb[(simulation_id, track_id)])
    return tad_hits_sample_dict

def convert_damage_file_in_tad_id_dict(damage_type, ion):
    damage_file_dict = defaultdict(lambda: np.array([]))
    damage_file_name = f"{PATH}/post_analysis/analysis_output/{ion}/Irradiation{IRRADIATION_ID}/{SIMULATED_DOSE}Gy/List_{damage_type}.txt"
    with open(damage_file_name, "r") as damage_file:
        [next(damage_file) for _ in range(0, 3)]
        for line in damage_file:
            if not line.strip():  #remove white spaces
                continue
            parts = line.split()
            run_id = parts[0]
            event_id = parts[1]
            tad_id = parts[3]
            damage_file_dict[(int(run_id), int(event_id))] = np.append(damage_file_dict[(int(run_id), int(event_id))], int(tad_id))
    return damage_file_dict

def write_tad_hit_file(damage_type, tad_hits, output_folder, nb_sample):
    name_file = output_folder + f"/tad_hits_{nb_sample}_{damage_type}.txt"
    with open(name_file, 'w') as file:
        for tad_id in tad_hits:
            file.write(f"{tad_id}\n")

def normalize_key(key):
    match = re.match(r"^[A-Za-z]+_\d+", key)
    return match.group(0) if match else None

def set_rng():
    np.random.seed(SEED)
    random.seed(SEED)

if __name__ == '__main__':
    set_rng()
    for ion in IONS:
        write_tad_hit_files_with_sampling(ion)
        print(f"{ion} is analysed\n")