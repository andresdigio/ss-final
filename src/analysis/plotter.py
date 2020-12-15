import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Current directory
curr_dir = os.getcwd()

# ../../data/
data_folder = os.path.abspath(os.path.join(curr_dir, os.pardir, os.pardir, 'data'))
data_files_paths = [os.path.join(data_folder, file_name) for file_name in os.listdir(data_folder)
                    if file_name.endswith('.data')]


def plot_count(df):
    plt.yticks(np.arange(0, 31, 2))
    plt.grid(b=True, which='both', axis='y', linestyle=':')
    plt.plot(df.t, df.n, label=df['N'][0])


def plot_total_energy(df):
    plt.plot(df.t, df.E, label=df['N'][0])


def plot_kinetic_energy(df):
    plt.plot(df.t, df.K, label=df['N'][0])


def plot_clockwise_particles(df):
    clockwise_percentage = df['clock']/df['n'] * 100
    plt.ylim(bottom=0, top=100)
    plt.plot(df.t, clockwise_percentage, label=df['N'][0])


def plot_collisions_per_window(df):
    collisions = df.loc[df['window'] != 0]
    t = df.t
    window = df.window
    plt.plot(t[window > 0], window[window > 0], linestyle='dotted', label=df['N'][0])


# Copy-pasted function, not sure what is going on but it sorts the legend labels in decreasing order
def get_handles_and_labels_for_sorted_legend():
    handles, labels = plt.gca().get_legend_handles_labels()
    return zip(*[(handles[i], labels[i]) for i in reversed(sorted(range(len(handles)),
                                                                  key=lambda k: list(map(int, labels))[k]))])


def plot_results(y_label, plot_function):
    plt.figure()
    plt.xlabel('Tiempo [s]')
    plt.ylabel(y_label)

    [plot_function(pd.read_csv(data_file)) for data_file in data_files_paths]

    handles, labels = get_handles_and_labels_for_sorted_legend()
    legend = plt.legend(handles, labels, loc='best', title='Partículas iniciales\nen el sistema')
    plt.setp(legend.get_title(), multialignment='center')
    plt.show()
    plt.close()


# plot_results('Energía total [J]', plot_total_energy)
# plot_results('Energía cinética [J]', plot_kinetic_energy)
# plot_results('% de partículas en sentido horario', plot_clockwise_particles)
plot_results('Colisiones por ventana', plot_collisions_per_window)
# plot_results('Cantidad de partículas', plot_count)
