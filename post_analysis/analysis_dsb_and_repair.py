import numpy as np
import os

# Open and read the file
file_name = "Output.dat"
ion_energy = "O_350_MeV"
output_folder = f"output/output{ion_energy}/"
nb_output_folder = 1000

def file_analysis():
    # write_tad_hit_files(output_folder)
    read_dsb_doses(output_folder)

def read_dsb_doses(output_folder):
    dsb = np.array([])
    dose = np.array([])
    for copy_nb in range(0, nb_output_folder):
        file = output_folder + f"Copy{copy_nb}/" + file_name
        dsb = np.append(dsb, extract_value_from_file(file, "DSB      [SB]", 2))
        dose = np.append(dose, extract_value_from_file(file, "Dose delivered in Cell Nucleus:", 5))
        print(f"Dsb analysis: {copy_nb + 1} / {nb_output_folder}")
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

def write_tad_hit_files(output_folder):
    tad_ion_folder = f"tad_hits/{ion_energy}"
    if not os.path.exists(tad_ion_folder):
        os.makedirs(tad_ion_folder)
    for copy_nb in range(0, nb_output_folder):
        output_folder_copy = output_folder + f"Copy{copy_nb}/"
        tad_hits_ssb = np.unique(get_tad_ids(output_folder_copy + "List_SSB.txt"))
        tad_hits_dsb = np.unique(get_tad_ids(output_folder_copy + "List_DSB.txt"))

        tad_hits_unique = np.unique(np.concatenate((tad_hits_ssb, tad_hits_dsb)))
        write_tad_hit_file(tad_ion_folder + f"/tad_hits_{copy_nb}_SSB.txt", tad_hits_ssb)
        write_tad_hit_file(tad_ion_folder + f"/tad_hits_{copy_nb}_DSB.txt", tad_hits_dsb)
        write_tad_hit_file(tad_ion_folder + f"/tad_hits_{copy_nb}_ALL.txt", tad_hits_unique)
        print(f"Tad hit analysis: {copy_nb + 1} / {nb_output_folder}")

def write_tad_hit_file(name_file, tad_hits):
    with open(name_file, 'w') as file:
        for tad_id in tad_hits:
            file.write(f"{tad_id}\n")

def get_tad_ids(filename):
    tad_ids = np.array([])
    with open(filename, 'r') as file:
        next(file)
        for line in file:
            parts = line.split()
            if len(parts) >= 2:
                tad_ids = np.append(tad_ids, int(parts[2]))
    return tad_ids


if __name__ == '__main__':
    file_analysis()
