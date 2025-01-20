import random
from collections import defaultdict
from scipy.stats import poisson
import numpy as np
import os.path

PATH = "/sps/gdrmi2b/levrague/dsb_and_repair"
IRRADIATION_ID = 1 #1 (HC nucleus) or 2 (decompacted nucleus)
NB_TRACKS_1_GY_DICT = {"O_350_MeV": 78}
ION_LET = {"O_350_MeV": 20.47} #keV/µm
DOSES = [1] #Gy
# DOSES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1] #Gy
NB_SIMULATIONS = 1000
NB_SIMULATION_SAMPLES = 1000
SEED = 42
ION = "O_350_MeV"

def nb_ions_poisson(dose):
    AREA = 16 * 16 #µm²
    UNIT_FACTOR = 6.24 #keV.kg.J^-1.µm^-3
    lambda_poisson = UNIT_FACTOR * dose * AREA / ION_LET[ION]
    #TO DO: use the real nb of ions per track to get x Gy (use dose output from Geant4 files)
    return poisson.rvs(lambda_poisson)

def write_tad_hit_files_with_sampling():
    damage_file_dict_ssb = convert_damage_file_in_tad_id_dict("SSB")
    damage_file_dict_dsb = convert_damage_file_in_tad_id_dict("DSB")
    for dose in DOSES:
        output_folder = f"{PATH}/post_analysis/tad_hits/{ION}/{dose}Gy"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for nb_sample in range(0, NB_SIMULATION_SAMPLES):
            nb_tracks_to_sample = nb_ions_poisson(dose)
            tad_hits_sample = compute_tad_hit_array_for_one_sample(damage_file_dict_ssb, damage_file_dict_dsb, nb_tracks_to_sample)
            write_tad_hit_file("SSB", np.unique(tad_hits_sample["SSB"]), output_folder, nb_sample)
            write_tad_hit_file("DSB", np.unique(tad_hits_sample["DSB"]), output_folder, nb_sample)
            write_tad_hit_file("ALL", np.unique(np.concatenate((tad_hits_sample["SSB"], tad_hits_sample["DSB"]))), output_folder, nb_sample)
            print(f"Sample n°{nb_sample + 1} / {NB_SIMULATION_SAMPLES}")
    return

def compute_tad_hit_array_for_one_sample(damage_file_dict_ssb, damage_file_dict_dsb, nb_tracks_to_sample):
    tad_hits_sample_dict = {"SSB": np.array([]), "DSB": np.array([]), "ALL": np.array([])}
    already_sampled_tracks = defaultdict(int)  # Simulation_id: track_id
    for nb_track in range(0, nb_tracks_to_sample):
        simulation_id = random.randint(0, NB_SIMULATIONS - 1)
        track_id = random.randint(0, NB_TRACKS_1_GY_DICT[ION])
        while (simulation_id, track_id) in already_sampled_tracks.items():
            simulation_id = random.randint(0, NB_SIMULATIONS - 1)
            track_id = random.randint(0, NB_TRACKS_1_GY_DICT[ION])
        already_sampled_tracks[simulation_id] = track_id
        tad_hits_sample_dict["SSB"] = np.append(tad_hits_sample_dict["SSB"], damage_file_dict_ssb[(simulation_id, track_id)])
        tad_hits_sample_dict["DSB"] = np.append(tad_hits_sample_dict["DSB"], damage_file_dict_dsb[(simulation_id, track_id)])
    return tad_hits_sample_dict

def convert_damage_file_in_tad_id_dict(damage_type):
    damage_file_dict = defaultdict(lambda: np.array([]))
    damage_file_name = f"{PATH}/post_analysis/analysis_output/{ION}/Irradiation{IRRADIATION_ID}/{DOSE}Gy/List_{damage_type}.txt"
    with open(damage_file_name, "r") as damage_file:
        next(damage_file)
        for line in damage_file:
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

def set_rng():
    np.random.seed(SEED)
    random.seed(SEED)

if __name__ == '__main__':
    set_rng()
    write_tad_hit_files_with_sampling()
