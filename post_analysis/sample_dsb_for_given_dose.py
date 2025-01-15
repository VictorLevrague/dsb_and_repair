import random
from collections import defaultdict
from scipy.stats import poisson
import numpy as np
import os.path

NB_TRACKS_1_GY_DICT = {"O_350_MeV": 78}
ION_LET = {"O_350_MeV": 20.47} #keV/µm
DOSES = [1] #Gy
# DOSES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1] #Gy
NB_SIMULATIONS = 10
NB_SIMULATION_SAMPLES = 10
SEED = 42

ion = "O_350_MeV"

def nb_ions_poisson(dose):
    AREA = 16 * 16 #µm²
    UNIT_FACTOR = 6.24 #keV.kg.J^-1.µm^-3
    lambda_poisson = UNIT_FACTOR * dose * AREA / ION_LET[ion]
    #TO DO: use the real nb of ions per track to get x Gy (use dose output from Geant4 files)
    return poisson.rvs(lambda_poisson)

def compute_nb_dsb():
    dsb_for_dose = {} #Dose[Gy]: DSB
    for dose in DOSES:
        nb_dsb = 0
        for nb_sample in range(0, NB_SIMULATION_SAMPLES):
            already_sampled_tracks = defaultdict(int) #Simulation_id: track_id
            nb_tracks_to_sample = nb_ions_poisson(dose)
            nb_dsb_one_Gy = 0
            for nb_track in range(0, nb_tracks_to_sample):
                simulation_id = random.randint(0, NB_SIMULATIONS - 1)
                track_id = random.randint(0, NB_TRACKS_1_GY_DICT[ion])
                while (simulation_id, track_id) in already_sampled_tracks.items():
                    simulation_id = random.randint(0, NB_SIMULATIONS - 1)
                    track_id = random.randint(0, NB_TRACKS_1_GY_DICT[ion])
                nb_dsb += compute_dsb_for_track(simulation_id, track_id)
                nb_dsb_one_Gy += compute_dsb_for_track(simulation_id, track_id)
                already_sampled_tracks[simulation_id] = track_id
            print(f"Sample n°{nb_sample + 1} / {NB_SIMULATION_SAMPLES}")
        dsb_for_dose[dose] = nb_dsb / NB_SIMULATION_SAMPLES
    return dsb_for_dose


def compute_dsb_for_track(copy_nb, track_id):
    file_name = f"output/output{ion}/Copy{copy_nb}/List_DSB.txt"
    nb_dsb = 0
    if not os.path.isfile(file_name):
        print("Non existing file")
        return nb_dsb
    with open(file_name, "r") as dsb_file:
        next(dsb_file)
        for line in dsb_file:
            parts = line.split()
            track_id_dsb_file = int(parts[0])
            if track_id_dsb_file == track_id:
                nb_dsb += 1
    return nb_dsb

def set_rng():
    np.random.seed(SEED)
    random.seed(SEED)

if __name__ == '__main__':
    set_rng()
    print(compute_nb_dsb())