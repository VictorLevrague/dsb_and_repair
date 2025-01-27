import random
import re
from collections import defaultdict
import pandas as pd
from scipy.stats import poisson
import numpy as np
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
#SAMPLE_DOSES = [1] #Gy
SAMPLE_DOSES = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1] #Gy
NB_SIMULATIONS = 1000
NB_SIMULATION_SAMPLES = 1000
SEED = 42
# "Cs_137_EC", "Cs_137_HC"
# IONS = ["H_250MeV_EC", "H_250MeV_HC",
#         "He_250MeV_EC", "He_250MeV_HC", "O_350MeV_EC", "O_350MeV_HC",
#         "O_55MeV_EC", "O_55MeV_HC", "Si_170MeV_EC", "Si_170MeV_HC",
#         "Ti_300MeV_EC", "Ti_300MeV_HC", "Fe_600MeV_EC", "Fe_600MeV_HC",
#         "Fe_450MeV_EC", "Fe_450MeV_HC", "Fe_300MeV_EC", "Fe_300MeV_HC"]
IONS = ["Fe_300MeV_EC", "Fe_300MeV_HC"]

def nb_ions_poisson(dose, ion):
    # AREA = 16 * 16 #µm²
    # UNIT_FACTOR = 6.24 #keV.kg.J^-1.µm^-3
    # lambda_poisson = UNIT_FACTOR * dose * AREA / ION_LET[ion]
    lambda_poisson = int(NB_TRACKS_1_GY_DICT[ion] * dose)
    return poisson.rvs(lambda_poisson)

def compute_nb_damage(damage_type, ion, damage_modifier):
    damage_for_dose = {} #Dose[Gy]: damage
    damage_for_dose_error = {}  # Dose[Gy]: damage standard deviation / square root N
    damage_file_dict = convert_damage_file_in_dict(damage_type, ion)
    for dose in SAMPLE_DOSES:
        damage_array = np.zeros(NB_SIMULATION_SAMPLES)
        for nb_sample in range(0, NB_SIMULATION_SAMPLES):
            already_sampled_tracks = defaultdict(int) #Simulation_id: track_id
            nb_tracks_to_sample = nb_ions_poisson(dose, normalize_key(ion))
            for _ in range(0, nb_tracks_to_sample):
                simulation_id = random.randint(0, NB_SIMULATIONS - 1)
                track_id = random.randint(0, NB_SIMULATED_TRACKS[normalize_key(ion)])
                while (simulation_id, track_id) in already_sampled_tracks.items():
                    simulation_id = random.randint(0, NB_SIMULATIONS - 1)
                    track_id = random.randint(0, NB_SIMULATED_TRACKS[normalize_key(ion)])
                damage_array[nb_sample] += damage_file_dict[(simulation_id, track_id)][damage_modifier]
                already_sampled_tracks[simulation_id] = track_id
            if (nb_sample+1)%100 ==0:
                print(f"Sample n°{nb_sample + 1} / {NB_SIMULATION_SAMPLES}")
        damage_for_dose[dose] = np.mean(damage_array)
        damage_for_dose_error[dose] = np.std(damage_array) / np.sqrt(NB_SIMULATION_SAMPLES)
    return damage_for_dose, damage_for_dose_error

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
            damage_file_dict[(int(run_id), int(event_id))]["NormalDamage"] += 1
            if int(parts[8]) > 1:
                damage_file_dict[(int(run_id), int(event_id))]["ComplexDamage"] += 1
            if int(parts[9]):
                damage_file_dict[(int(run_id), int(event_id))]["DirectDamage"] += 1
            if int(parts[10]):
                damage_file_dict[(int(run_id), int(event_id))]["IndirectDamage"] += 1
    return damage_file_dict

def normalize_key(key):
    match = re.match(r"^[A-Za-z]+_\d+", key)
    return match.group(0) if match else None

def set_rng():
    np.random.seed(SEED)
    random.seed(SEED)

if __name__ == '__main__':
    set_rng()
    output_dict = defaultdict(lambda: np.array([]))
    for ion in IONS:
        ssb, ssb_error = compute_nb_damage("SSB", ion, "NormalDamage")
        ssb_direct, ssb_direct_error = compute_nb_damage("SSB", ion, "DirectDamage")
        ssb_indirect, ssb_indirect_error = compute_nb_damage("SSB", ion, "IndirectDamage")
        dsb, dsb_error = compute_nb_damage("DSB", ion, "NormalDamage")
        dsb_direct, dsb_direct_error = compute_nb_damage("DSB", ion, "DirectDamage")
        dsb_indirect, dsb_indirect_error = compute_nb_damage("DSB", ion, "IndirectDamage")
        dsb_complex, dsb_complex_error = compute_nb_damage("DSB", ion, "ComplexDamage")
        for dose in ssb.keys():
            output_dict["Ion"] = np.append(output_dict["Ion"], ion)
            output_dict["LET"] = np.append(output_dict["LET"], ION_LET[normalize_key(ion)])
            output_dict["Dose"] = np.append(output_dict["Dose"], dose)
            output_dict["SSB"] = np.append(output_dict["SSB"], ssb[dose])
            output_dict["SSB_error"] = np.append(output_dict["SSB_error"], ssb_error[dose])
            output_dict["SSB_direct"] = np.append(output_dict["SSB_direct"], ssb_direct[dose])
            output_dict["SSB_direct_error"] = np.append(output_dict["SSB_direct_error"], ssb_direct_error[dose])
            output_dict["SSB_indirect"] = np.append(output_dict["SSB_indirect"], ssb_indirect[dose])
            output_dict["SSB_indirect_error"] = np.append(output_dict["SSB_indirect_error"], ssb_indirect_error[dose])
            output_dict["DSB"] = np.append(output_dict["DSB"], dsb[dose])
            output_dict["DSB_error"] = np.append(output_dict["DSB_error"], dsb_error[dose])
            output_dict["DSB_direct"] = np.append(output_dict["DSB_direct"], dsb_direct[dose])
            output_dict["DSB_direct_error"] = np.append(output_dict["DSB_direct_error"], dsb_direct_error[dose])
            output_dict["DSB_indirect"] = np.append(output_dict["DSB_indirect"], dsb_indirect[dose])
            output_dict["DSB_indirect_error"] = np.append(output_dict["DSB_indirect_error"], dsb_indirect_error[dose])
            output_dict["DSB_complex"] = np.append(output_dict["DSB_complex"], dsb_complex[dose])
            output_dict["DSB_complex_error"] = np.append(output_dict["DSB_complex_error"], dsb_complex_error[dose])
            print(f"{dose} is analysed")
        print(f"{ion} is analysed\n")
    output_df = pd.DataFrame(output_dict)
    output_df.to_csv("damage_sample_analysis.csv", index=False)
