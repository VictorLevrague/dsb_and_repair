import numpy as np
import os

ION_ENERGY = "O_350_MeV"
NB_SIMULATIONS = 1000

def merge_ssb_dsb_files():
    analysis_output_folder = f"analysis_output/{ION_ENERGY}"
    if not os.path.exists(analysis_output_folder):
        os.makedirs(analysis_output_folder)
    write_output_file("SSB", analysis_output_folder)
    write_output_file("DSB", analysis_output_folder)

def write_output_file(damage_type, output_folder):
    with open(f"{output_folder}/List_{damage_type}.txt", "w") as output_file:
        output_file.write("0:run\t1:event\t2:chromosome\t3:domain\t4:genomic_position(bp)\t5:x(nm)\t6:y(nm)\t7:z(nm)\t8:complexity\t9:direct\t10:indirect\n")
        for copy_nb in range(0, NB_SIMULATIONS):
            damage_file_name = f"output/output{ION_ENERGY}/Copy{copy_nb}/List_{damage_type}.txt"
            with open(damage_file_name) as damage_file:
                next(damage_file)
                for line in damage_file:
                    output_file.write(f"{copy_nb}\t" + line)
            print(f"Merged file n°{copy_nb+1} / {NB_SIMULATIONS} for {damage_type}")

if __name__ == '__main__':
    merge_ssb_dsb_files()