import numpy as np
import os

# Open and read the file
FILE_NAME = "Output.dat"
ION_ENERGY = "O_350_MeV"
NB_OUTPUT_FOLDER = 1000
PATH = "/sps/gdrmi2b/levrague/dsb_and_repair"
IRRADIATION_ID = 1 #1 (HC nucleus) or 2 (decompacted nucleus)

def file_analysis():
    output_folder = f"{PATH}/output/output{ION_ENERGY}/outputIrradiation{IRRADIATION_ID}/"
    read_dsb_doses(output_folder)

def read_dsb_doses(output_folder):
    dsb = np.array([])
    dose = np.array([])
    for copy_nb in range(0, NB_OUTPUT_FOLDER):
        file = output_folder + f"Copy{copy_nb}/" + FILE_NAME
        dsb = np.append(dsb, extract_value_from_file(file, "DSB      [SB]", 2))
        dose = np.append(dose, extract_value_from_file(file, "Dose delivered in Cell Nucleus:", 5))
        print(f"Dsb analysis: {copy_nb + 1} / {NB_OUTPUT_FOLDER}")
    print("average DSB: ", np.mean(dsb), " +- ", np.std(dsb))
    print("average Dose: ", np.mean(dose), " +- ", np.std(dose))
    print("average DSB/Gy: ", np.mean(dsb) / np.mean(dose))
    evaluate_standard_error(dsb, dose)
    return (dsb, dose)

def extract_value_from_file(file_path, value_type, string_index):
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip().startswith(value_type):
                parts = line.split()
                if len(parts) > 1:
                    value = float(parts[string_index])
                    return value
    print(f"Value not found in file: {file_path}")
    return 0

def bootstrap_ratios(breaks, doses, n_bootstrap):
    bootstrap_ratios = []
    for _ in range(n_bootstrap):
        indices = np.random.choice(len(breaks), size=len(breaks), replace=True)
        resampled_breaks = breaks[indices]
        resampled_doses = doses[indices]
        bootstrap_ratios.append(np.sum(resampled_breaks) / np.sum(resampled_doses))
    return np.array(bootstrap_ratios)

def evaluate_standard_error(breaks, doses):
    n_bootstrap_values = [10000]
    for n in n_bootstrap_values:
        bootstrap_results = bootstrap_ratios(breaks, doses, n)
        error_estimate = np.std(bootstrap_results)
        print(f"Bootstrap samples: {n}, Standard Deviation estimate: {error_estimate}")
        print(f"Bootstrap samples: {n}, Standard Error estimate: {error_estimate/np.sqrt(len(breaks))}")
        print(f"Bootstrap samples: {n}, 2* Standard Error estimate: {2*error_estimate / np.sqrt(len(breaks))}")


if __name__ == '__main__':
    file_analysis()