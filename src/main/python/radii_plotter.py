import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from plotting_utils import get_handles_and_labels_for_sorted_legend, get_data_folder, plot_mean_and_error, JAR_PATH


def plot_radii_with_mean_and_error(radii):
    plt.figure()
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Cantidad de partículas orbitando')

    plot_mean_and_error(radii)

    plt.yticks(np.arange(0, 21, 2))
    handles, labels = get_handles_and_labels_for_sorted_legend()
    legend = plt.legend(handles, labels, loc='best', title='Radio de partículas [m]')
    plt.setp(legend.get_title(), multialignment='center')
    plt.show()
    plt.close()


def run_simulation(T, dt, o, radius):
    cmd = 'java -jar ' + JAR_PATH + ' -data -T {:g} -dt {:g} -o {:g} -r {:g}'.format(T, dt, o, radius)
    print(cmd)
    os.system(cmd)


def run_multiple_simulations(runs, radii):
    for radius in radii:
        for run in runs.tolist():
            run_simulation(2500, 5e-5, 0.7, radius)
            move_data_file_to_own_folder(radius, run)


def move_data_file_to_own_folder(radius, run):
    data_folder = get_data_folder()
    new_data_file_name = [file_name for file_name in os.listdir(data_folder) if file_name.endswith('.data')][0]
    old_file_path = data_folder + '/' + new_data_file_name
    new_folder = data_folder + '/r' + str(radius)

    if not os.path.isdir(new_folder):
        os.mkdir(new_folder)

    new_file_path = new_folder + '/' + new_data_file_name.split('.data')[0] + '(' + str(run) + ')' + '.data'
    os.replace(old_file_path, new_file_path)


radii_arg = [1, 2, 3]
runs_arg = np.arange(5)
run_multiple_simulations(runs_arg, radii_arg)
plot_radii_with_mean_and_error(radii_arg)


# Not sure if these might needed

# def get_radius(name):
#     return int(((name.split(',')[-1]).split('.data')[0]).split('R')[-1])
#
#
# def plot_results(y_label, plot_function):
#     plt.figure()
#     plt.xlabel('Tiempo [s]')
#     plt.ylabel(y_label)
#     [plot_function(pd.read_csv(data_file), get_radius((data_file.split('/'))[-1])) for data_file in data_files_paths]
#     handles, labels = get_handles_and_labels_for_sorted_legend()
#     legend = plt.legend(handles, labels, loc='best', title='Radios de las partículas [m]')
#     plt.setp(legend.get_title(), multialignment='center')
#     plt.show()
#     plt.close()
