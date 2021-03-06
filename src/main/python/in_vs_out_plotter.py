import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from plotting_utils import get_data_files, get_handles_and_labels_for_sorted_legend, get_data_folder, average, JAR_PATH


def plot_in_and_out_with_mean_and_error(N_particles):
    plt.figure()
    plt.xlabel('% inicial de partículas en sentido horario')
    plt.ylabel('% final de partículas en sentido horario')

    # Plotting logic
    for N in N_particles:
        N_folder = get_data_folder('in_vs_out/n' + str(N))
        in_o = []
        out_o = []

        for percentage in os.listdir(N_folder):
            percentage_folder = N_folder + '/' + percentage
            data_frames = [pd.read_csv(os.path.join(percentage_folder, data_file)) for data_file in os.listdir(percentage_folder)]
            in_o.append(100 * average([df.o.iat[0] for df in data_frames]))
            out_o.append(100 * average([df.clock.iat[-1]/df.n.iat[-1] for df in data_frames]))

        plt.plot(in_o, out_o, marker='.', mfc='c', label=N)

    plt.xticks(np.arange(0, 101, 10))
    plt.yticks(np.arange(0, 101, 10))
    handles, labels = get_handles_and_labels_for_sorted_legend()
    legend = plt.legend(handles, labels, loc='best', title='Cantidad inicial de partículas')
    plt.setp(legend.get_title(), multialignment='center')
    plt.show()
    plt.close()


def run_simulation(T, dt, o, N):
    cmd = 'java -jar ' + JAR_PATH + ' -data -T {:g} -dt {:g} -o {:g} -N {:g}'.format(T, dt, o, N)
    print(cmd)
    os.system(cmd)


def run_multiple_simulations(runs, N_particles):
    for N in N_particles:
        for run in runs:
            run_simulation(2500, 5e-5, 0.7, N)
            move_data_file_to_own_folder(N, run)


def move_data_file_to_own_folder(N, run):
    data_folder = get_data_folder()
    new_data_file_name = [file_name for file_name in os.listdir(data_folder) if file_name.endswith('.data')][0]
    old_file_path = data_folder + '/' + new_data_file_name
    new_folder = data_folder + '/N' + str(N)

    if not os.path.isdir(new_folder):
        os.mkdir(new_folder)

    new_file_path = new_folder + '/' + new_data_file_name.split('.data')[0] + '(' + str(run) + ')' + '.data'
    os.replace(old_file_path, new_file_path)


N_particles_arg = [10, 15, 20]
runs_arg = np.arange(5)
# run_multiple_simulations(runs_arg, N_particles_arg)
plot_in_and_out_with_mean_and_error(N_particles_arg)
# plot_in_and_out_with_mean_and_error([10])
# plot_in_and_out_with_mean_and_error([15])
# plot_in_and_out_with_mean_and_error([20])
