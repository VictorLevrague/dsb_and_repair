from collections import defaultdict
import numpy as np
import os

PATH = "/sps/gdrmi2b/levrague/dsb_and_repair"
ION = "O_350_MeV"
DOSE = 1 #Gy
NB_OUTPUT_FOLDER = 1000

def file_analysis():
    output_folder = f"{PATH}/post_analysis/tad_hits/{ION}/{DOSE}Gy"
    write_tad_hit_files(output_folder)

def write_tad_hit_files(output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    damage_file_dict_ssb = convert_damage_file_in_tad_id_dict("SSB")
    damage_file_dict_dsb = convert_damage_file_in_tad_id_dict("DSB")
    for copy_nb in range(0, NB_OUTPUT_FOLDER):
        tad_hits_ssb = np.sort(np.array([tad_id for (run_id, tad_id) in damage_file_dict_ssb.keys() if run_id == copy_nb ]))
        tad_hits_dsb = np.sort(np.array([tad_id for (run_id, tad_id) in damage_file_dict_dsb.keys() if run_id == copy_nb ]))
        tad_hits_unique = np.unique(np.concatenate((tad_hits_ssb, tad_hits_dsb)))
        write_tad_hit_file(output_folder + f"/tad_hits_{copy_nb}_SSB.txt", tad_hits_ssb)
        write_tad_hit_file(output_folder + f"/tad_hits_{copy_nb}_DSB.txt", tad_hits_dsb)
        write_tad_hit_file(output_folder + f"/tad_hits_{copy_nb}_ALL.txt", tad_hits_unique)
        print(f"Tad hit analysis: {copy_nb + 1} / {NB_OUTPUT_FOLDER}")

def write_tad_hit_file(name_file, tad_hits):
    with open(name_file, 'w') as file:
        for tad_id in tad_hits:
            file.write(f"{tad_id}\n")

def convert_damage_file_in_tad_id_dict(damage_type):
    damage_file_dict = defaultdict(int)
    damage_file_name = f"{PATH}/post_analysis/analysis_output/{ION}/Irradiation1/{DOSE}Gy/List_{damage_type}.txt"
    with open(damage_file_name, "r") as damage_file:
        next(damage_file)
        for line in damage_file:
            parts = line.split()
            run_id = parts[0]
            tad_id = parts[3]
            damage_file_dict[(int(run_id), int(tad_id))] = 1
    return damage_file_dict

if __name__ == '__main__':
    file_analysis()
