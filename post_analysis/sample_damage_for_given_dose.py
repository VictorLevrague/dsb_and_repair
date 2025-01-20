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

ion = "O_350_MeV"
IRRADIATION_ID = 1 #1 (HC nucleus) or 2 (decompacted nucleus)
DOSE = 1 #Gy

def nb_ions_poisson(dose):
    AREA = 16 * 16 #µm²
    UNIT_FACTOR = 6.24 #keV.kg.J^-1.µm^-3
    lambda_poisson = UNIT_FACTOR * dose * AREA / ION_LET[ion]
    #TO DO: use the real nb of ions per track to get x Gy (use dose output from Geant4 files)
    return poisson.rvs(lambda_poisson)

def compute_nb_damage(damage_type):
    damage_for_dose = {} #Dose[Gy]: damage
    damage_file_dict = convert_damage_file_in_dict(damage_type)
    for dose in DOSES:
        nb_damage = 0
        for nb_sample in range(0, NB_SIMULATION_SAMPLES):
            already_sampled_tracks = defaultdict(int) #Simulation_id: track_id
            nb_tracks_to_sample = nb_ions_poisson(dose)
            for nb_track in range(0, nb_tracks_to_sample):
                simulation_id = random.randint(0, NB_SIMULATIONS - 1)
                track_id = random.randint(0, NB_TRACKS_1_GY_DICT[ion])
                while (simulation_id, track_id) in already_sampled_tracks.items():
                    simulation_id = random.randint(0, NB_SIMULATIONS - 1)
                    track_id = random.randint(0, NB_TRACKS_1_GY_DICT[ion])
                nb_damage += damage_file_dict[(simulation_id, track_id)]
                already_sampled_tracks[simulation_id] = track_id
            print(f"Sample n°{nb_sample + 1} / {NB_SIMULATION_SAMPLES}")
        damage_for_dose[dose] = nb_damage / NB_SIMULATION_SAMPLES
    return damage_for_dose

def convert_damage_file_in_dict(damage_type):
    damage_file_dict = defaultdict(int)
    damage_file_name = f"{PATH}/post_analysis/analysis_output/{ion}/Irradiation{IRRADIATION_ID}/{DOSE}Gy/List_{damage_type}.txt"
    with open(damage_file_name, "r") as damage_file:
        next(damage_file)
        for line in damage_file:
            parts = line.split()
            run_id = parts[0]
            event_id = parts[1]
            damage_file_dict[(int(run_id), int(event_id))] += 1
    return damage_file_dict

def set_rng():
    np.random.seed(SEED)
    random.seed(SEED)

if __name__ == '__main__':
    set_rng()
    print(compute_nb_damage("DSB"))
