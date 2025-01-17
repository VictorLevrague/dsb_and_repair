import numpy as np
import os
import time
from collections import defaultdict

PATH = "/sps/gdrmi2b/levrague/dsb_and_repair"
ION_ENERGY = "O_350_MeV"
DOSE = 1 #Gy
NB_OUTPUT_FOLDER = 10

def write_geometry_file(output_filename, tad_hits, header, hc_geom_dict, ec_geom_dict):
    hit_tads = read_hit_tads(tad_hits)
    with open(output_filename, 'w') as file:
        file.writelines(header)
        hit_tads_set = set(hit_tads) #slightly faster search time
        nb_tads = 15282
        for tad_id in range(0, nb_tads):
            # print(f"Writing TAD {tad_id+1} / {nb_tads}")
            if tad_id in hit_tads_set:
                file.writelines(ec_geom_dict[tad_id])
            else:
                file.writelines(hc_geom_dict[tad_id])

def read_geometry_files(filename):
    with open(filename, 'r') as file:
        header = [file.readline() for _ in range(11)]
        tad_data = file.readlines()
    print(f"Geometry file {filename} read")
    return header, tad_data

def group_by_tad_id(tad_data):
    tad_dict = defaultdict(list)
    for line in tad_data:
        parts = line.split()
        if len(parts) >= 4:
            tad_id = int(parts[3])
            tad_dict[tad_id].append(line)
    print("Geometry grouped by TAD")
    return tad_dict

def read_hit_tads(namefile):
    hit_tads = np.array([])
    with open(namefile, 'r') as file:
        for line in file:
            hit_tads = np.append(hit_tads, int(float(line.split()[0])))
    return hit_tads

if __name__ == "__main__":
    start_time1 = time.time()
    header, voxels_tad_hc = read_geometry_files("geometrie_fibroblast_hc_floriane.fab2g4dna")
    voxels_tad_ec = read_geometry_files("geometrie_fibroblast_ec_floriane.fab2g4dna")[1]
    voxels_tad_hc_dict = group_by_tad_id(voxels_tad_hc)
    voxels_tad_ec_dict = group_by_tad_id(voxels_tad_ec)
    ion_folder = f"output/output{ION_ENERGY}/{DOSE}Gy"
    if not os.path.exists(ion_folder):
        os.makedirs(ion_folder)
    start_time2 = time.time()
    print("time before loop: ", start_time2 - start_time1, " s")
    for copy_nb in range(0, NB_OUTPUT_FOLDER):
        write_geometry_file(ion_folder + f"/irradiated_geometry_file_{copy_nb}.fab2g4dna", f"{PATH}/post_analysis/tad_hits/{ION_ENERGY}/{DOSE}Gy/tad_hits_{copy_nb}_ALL.txt", header, voxels_tad_hc_dict, voxels_tad_ec_dict)
        print(f"Geometry n°{copy_nb+1} created")
    print("time per output: ", (time.time()-start_time2) / NB_OUTPUT_FOLDER, " s")