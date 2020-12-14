import os
import pandas as pd
import matplotlib.pyplot as plt

# Current directory
curr_dir = os.getcwd()

# ../../data/
data_folder = os.path.abspath(os.path.join(curr_dir, os.pardir, os.pardir, 'data'))
data_files_paths = [os.path.join(data_folder, file_name) for file_name in os.listdir(data_folder)
                    if file_name.endswith('.data')]

plt.figure()
plt.xlabel('Tiempo [s]')
plt.ylabel('Energía [J]')

for data_file in data_files_paths:
    df = pd.read_csv(data_file)
    energy_time = df[['t', 'E']]
    plt.plot(energy_time.t, energy_time.E)

plt.show()


plt.figure()
plt.xlabel('Tiempo [s]')
plt.ylabel('Partículas en sentido horario')

for data_file in data_files_paths:
    df = pd.read_csv(data_file)
    clockwise_particles = df[['t', 'clock']]
    plt.plot(clockwise_particles.t, clockwise_particles.clock)

plt.show()

