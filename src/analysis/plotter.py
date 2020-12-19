import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Current directory
curr_dir = os.getcwd()


def get_data_folder(sub_folder=''):
    return os.path.abspath(os.path.join(curr_dir, os.pardir, os.pardir, 'data/' + sub_folder))


def get_data_files(sub_folder):
    # ../../data/ + sub_folder
    data_folder = get_data_folder(sub_folder)
    return [os.path.join(data_folder, file_name) for file_name in os.listdir(data_folder)
            if file_name.endswith('.data')]


def plot_survivors(N=20):
    df = pd.read_csv(get_data_folder() + '/survivors.csv')
    df = df[df['N'] == N]
    df['o'] = df.apply(lambda x: x.o*100, axis=1)
    df['s_rate'] = df.apply(lambda x: x.s*100/x.N, axis=1)
    plt.figure()
    gp = df.groupby(['o'])
    means = gp.mean()
    error = gp.std()

    plt.xlabel('Distribucion de orientacion inicial en sentido horario [%]')
    plt.ylabel('Porcentaje de supervivencia en equilibrio [%]')
    plt.errorbar(means.index, means['s_rate'], yerr=error['s_rate'], marker='o', fmt='o', mfc='c', capthick=10, ms=7)

    plt.savefig('C:/Users/Andres/ss-final/graphs/{:d}_survival.png'.format(N))
    plt.show()


def plot_energy_dt(df):
    plt.plot(df.t, df.E, label=df['dt'][0])


def plot_count(df):
    plt.yticks(np.arange(0, 51, 2))
    plt.grid(b=True, which='both', axis='y', linestyle=':')
    plt.plot(df.t, df.n, label=df['N'][0])


def plot_total_energy(df):
    plt.plot(df.t, df.E, label=df['N'][0])


def plot_kinetic_energy(df):
    plt.plot(df.t, df.K, label=df['N'][0])


def plot_clockwise_particles(df):
    plt.ylim(bottom=0, top=105)
    clockwise_percentage = df['clock']/df['n'] * 100
    label = '{:d}'.format(int(df['o'].mode()[0]*100))
    plt.plot(df.t, clockwise_percentage, label=label, linestyle='dotted')


def plot_collisions(df):
    plt.scatter(df.t, df.collisions, linestyle='-', label=df['N'][0], marker='.')


# Copy-pasted function, not sure what is going on but it sorts the legend labels in decreasing order
def get_handles_and_labels_for_sorted_legend():
    handles, labels = plt.gca().get_legend_handles_labels()
    return zip(*[(handles[i], labels[i]) for i in reversed(sorted(range(len(handles)),
                                                                  key=lambda k: list(map(int, labels))[k]))])


def plot_results(y_label, plot_function):
    plt.figure()
    plt.xlabel('Tiempo [s]')
    plt.ylabel(y_label)
    [plot_function(pd.read_csv(data_file)) for data_file in get_data_files()]
    handles, labels = get_handles_and_labels_for_sorted_legend()
    legend = plt.legend(handles, labels, loc='best', title='Proporción inicial\nen el sistema')
    plt.setp(legend.get_title(), multialignment='center')
    # plt.savefig('C:/Users/Andres/ss-final/graphs/stationary.png')
    plt.show()
    plt.close()


def calc_df_sum(dataframes, col_name):
    for df in dataframes[1:]:
        dataframes[0][col_name] += df[col_name]
    return dataframes[0]


def plot_count_with_mean_and_error():
    plt.figure()
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Cantidad de partículas')

    for percentage in np.arange(50, 51, 20):
        data_frames = [pd.read_csv(data_file) for data_file in get_data_files(str(percentage))]
        gp = pd.concat(data_frames).groupby(['t'])
        means = gp.mean()
        error = gp.std()
        plt.yticks(np.arange(0, 51, 2))
        plt.grid(b=True, which='both', axis='y', linestyle=':')
        plt.errorbar(means.index, means['n'], yerr=error['n'], marker='.', fmt='o', mfc='c', capthick=2, ms=7, label=data_frames[0]['N'][0])
    handles, labels = get_handles_and_labels_for_sorted_legend()
    legend = plt.legend(handles, labels, loc='best', title='% de partículas en sentido horario')
    plt.setp(legend.get_title(), multialignment='center')
    plt.show()
    plt.close()


# plot_results('Energía total [J]', plot_total_energy)
# plot_results('Energía cinética [J]', plot_kinetic_energy)
# plot_results('% de partículas en sentido horario', plot_clockwise_particles)
# plot_results('Colisiones', plot_collisions)
# plot_results('Cantidad de partículas', plot_count)
plot_count_with_mean_and_error()
# plot_survivors(N=50)
