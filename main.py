################################################################################################***** ***
# This program calculates Kubelka-Munk Function and data for Tauc plots from
# reflectance UV-Vis reflectance spectra, measured on lambda 750 in 1.11 laboratory,
# (through the door and to the left), in ACMiN, Kraków, Poland, Europe, Earth, Solar System, Milky Way.
# Results can be saved to txt file for further processing and visualisation. Bandgap can be estimated
# through line fitting of chosen tauc plot.
################################################################################################***** ***

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

# Necessary variables so that program doesn't crash
list_of_csvfiles = []
band_gap = ""
collection_of_bandgaps = {}
tauc = []

print("***** ***\nThis program calculates Kubelka-Munk Function, Tauc plots and Bandgaps from \n "
      "reflectance UV-Vis reflectance spectra, measured on lambda 750 in 1.11 laboratory, \n "
      "(through the door and to the left), in ACMiN, Kraków, Poland, Europe, Earth, Solar System, Milky Way.\n"
      "Results can be saved to txt file for further processing and visualisation. Bandgap can be estimated\n"
      "through line fitting of chosen tauc plot.\n"
      "***** ***\n")

# Working directory
print("Enter absolute path of directory with uv-vis data in X:\Y\Z format (single backlash)")
current_working_directory = input()

try:
    directory = os.chdir(current_working_directory)
except:
    print("No such directory exists")
    input("Press enter to quit")
    quit()

var1, var2 = [int(x) for x in
              input("\nEnter wavelenght range separated with space: minimum maximum (e.g. 200 1500)\n").split()]

saving_txt = input("\nDo you want to save results as txt file? (y/n) \n")

plotting = input("\nDo you want to plot results? (y/n) \n")
if plotting == "y":
    band_gap = input("\nDo you want to determine band gap? (y/n) \n")
    if band_gap == "y":
        transition_type = input(
            '\nWhat type of transition do you want to determine? 1. Direct allowed or 2. Indirect allowed? (press 1 or '
            '2) \n')

data_range = np.arange(var1, var2 + 1)
hv = 1240 / data_range

# Defining all necesary function for data transformation, plotting and band gap estimation
def data_transformation(x):
    # Reading data file
    global tauc1, tauc2
    data_raw = np.genfromtxt(x, delimiter=',', skip_header=2)

    # data preparation
    data = np.flip(data_raw[:, 1], 0)

    # data normalisation
    max_arg = np.amax(data)
    norm_data = data / max_arg

    # Kubelk-Munk Function
    Fkm = ((1 - norm_data) ** 2) / (2 * norm_data)

    # Preparation for Tauc plots
    try:
        tauc1 = (Fkm * hv) ** 2  # direct allowed
        tauc2 = (Fkm * hv) ** (1 / 2)  # indirect allowed
    except:
        print("\nWrong wavelenght range.")
        input("Press enter to quit")
        quit()
    return data_raw, data, norm_data, Fkm, tauc1, tauc2

def plot(norm_data, Fkm, tauc1, tauc2):
    figure, axis = plt.subplots(1, 4, squeeze=False)
    axis[0, 0].plot(data_range, norm_data)
    axis[0, 0].set_title("Normalised reflectance")
    axis[0, 1].plot(data_range, Fkm)
    axis[0, 1].set_title("Kubelka-Munk Function")
    axis[0, 2].plot(hv, tauc1)
    axis[0, 2].set_title("Tauc plot, direct allowed")
    axis[0, 3].plot(hv, tauc2)
    axis[0, 3].set_title("Tauc plot, indirect allowed")
    plt.show()

def get_bandgap(hv, tauc):
    hv_min, hv_max = [float(x) for x in input(
        "\nEnter energy range for line fitting in tauc plot: minimum maximum (e.g. 3.2 3.8).\n").split()]

    # Energy scale to wavelenght for indexing. Also values correction (offset from starting wavelenght)
    fit_min = int(1240.0 / hv_max) - var1
    fit_max = int(1240.0 / hv_min) - var1

    # Variable correction
    tauc = tauc[0]

    # Line fitting and band gap value extraction from X-axis intercept
    x = hv[fit_min:fit_max]
    y = tauc[fit_min:fit_max]
    linear_model = np.polyfit(x, y, 1)
    linear_model_fn = np.poly1d(linear_model)
    band_gap_value = (-linear_model[1]) / (linear_model[0])

    # Plotting with fitted line
    plt.xlabel("Energy / eV")
    plt.plot(hv, np.zeros(len(hv)), color="black")
    plt.plot(hv, linear_model_fn(hv), color="red")
    plt.plot(hv, tauc)
    plt.show()
    return band_gap_value

# Iteration over all csv files in working directory
for file in os.listdir(directory):
    if file.endswith(".csv"):
        file_name = str(file)
        list_of_csvfiles.append(file_name)

        data_raw, data, norm_data, Fkm, tauc1, tauc2 = data_transformation(file_name)
        # Plotting
        if plotting == "y":
            plot(norm_data, Fkm, tauc1, tauc2)
        else:
            pass

        # Bandgap calculation and plotting
        if band_gap == "y" and transition_type == "1":
            tauc.append(tauc1)
        else:
            tauc.append(tauc2)

        if band_gap == "y":
            band_gap_value = get_bandgap(hv, tauc)
            band_gap_value_listed = [band_gap_value]
            collection_of_bandgaps[file_name] = band_gap_value_listed
            band_gap_value_listed = []
        else:
            pass

        if saving_txt == "y":
            # Preparation for array construction
            mega_array = np.stack((data_range, hv, data, norm_data, Fkm, tauc1, tauc2))
            mega_array_T = np.transpose(mega_array)
            column_names = ["wavelenght", "energy", "raw data", "normalised data", "Fkm", "tauc^2 (direct allowed)",
                            "tauc^(1/2) (indirect allowed)"]
            df = pd.DataFrame(data=mega_array_T, columns=column_names)

            # Save file as txt
            df.to_csv(file_name.replace(".csv", '.txt'), sep='\t', index=False)
        else:
            pass

        tauc = []

# Raise error if wrong directory was chosen
if not list_of_csvfiles:
    input("\nNo csv files in chosen directory. Run program again with correct path")

saving_bandaps = pd.DataFrame.from_dict(collection_of_bandgaps)
saving_bandaps_T = saving_bandaps.T
saving_bandaps_T.to_csv("collection_of_bandgaps.txt", sep='\t')

input("\nEverything is done! Press Enter to quit and have nice day.")

