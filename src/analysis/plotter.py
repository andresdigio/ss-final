import os
import pandas as pd
import matplotlib.pyplot as plt

# Current directory
curr_dir = os.getcwd()

# ../../data/
data_folder = os.path.abspath(os.path.join(curr_dir, os.pardir, os.pardir, 'data'))
data_files_paths = [os.path.join(data_folder, file_name) for file_name in os.listdir(data_folder)
                    if file_name.endswith('.data')]

def plot_energy():
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Energía [J]')
    for data_file in data_files_paths:
        df = pd.read_csv(data_file)
        n = data_file.split(',')[-1].split('.')[0]
        energy_time = df[['t', 'E']]
        plt.plot(energy_time.t, energy_time.E, label=n)

    plt.legend(loc='best')

    plt.show()

def plot_orientation():
    plt.figure()
    plt.xlabel('Tiempo [s]')
    plt.ylabel('% de partículas en sentido horario')
    plt.ylim((0, 100))
    for data_file in data_files_paths:
        df = pd.read_csv(data_file)
        prcnt = df.apply(lambda row: 100*row['clock']/row['n'], axis=1)
        clockwise_particles = df[['t', 'clock']]
        plt.plot(clockwise_particles.t, prcnt)

    plt.show()

def plot_count():
    plt.figure()
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Cantidad de particulas')
    for data_file in data_files_paths:
        df = pd.read_csv(data_file)
        plt.plot(df.t, df.n)

    plt.show()

plot_orientation()
plot_energy()
plot_count()

