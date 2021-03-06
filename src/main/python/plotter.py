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
    df['cp'] = df['clock']/df['n'] * 100
    gp = df.groupby(['t'])
    means = gp.mean()
    error = gp.std()
    label = '{:d}'.format(int(df['o'].mode()[0]*100))
    plt.errorbar(means.index, means.cp, yerr=error.cp, label=label, fmt='o', marker='o', capthick=10, ms=3)

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
    legend = plt.legend(handles, labels, loc='lower left', title='Proporción inicial\nen el sistema')
    plt.setp(legend.get_title(), multialignment='center')
    plt.show()
    plt.close()


def calc_df_sum(dataframes, col_name):
    for df in dataframes[1:]:
        dataframes[0][col_name] += df[col_name]
    return dataframes[0]


def plot_count_with_mean_and_error(orientation_percentages):
    plt.figure()
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Cantidad de particulas')

    for percentage in orientation_percentages:
        data_frames = [pd.read_csv(data_file) for data_file in get_data_files(str(percentage))]
        gp = pd.concat(data_frames).groupby(['t'])
        means = gp.mean()
        error = gp.std()
        plt.yticks(np.arange(0, 51, 2))
        plt.grid(b=True, which='both', axis='y', linestyle=':')
        plt.errorbar(means.index, means['n'], yerr=error['n'], marker='.', fmt='o', mfc='c', capthick=2, ms=7, label=percentage)
    handles, labels = get_handles_and_labels_for_sorted_legend()
    legend = plt.legend(handles, labels, loc='best', title='% de partículas')
    plt.setp(legend.get_title(), multialignment='center')
    plt.show()
    plt.close()

def plot_time_with_mean_and_error(orientation_percentages):
    plt.figure()
    plt.ylabel('Tiempo [s]')
    plt.xlabel('Distribucion de orientacion inicial en sentido horario [%]')

    df = pd.DataFrame(orientation_percentages)

    for percentage in orientation_percentages:
        data_frames = [pd.read_csv(data_file) for data_file in get_data_files(str(percentage))]
        times = [df.iloc[-1].t for df in data_frames]
        tmpdf = pd.DataFrame(times, columns=['t'])
        tmpdf['o'] = percentage

        df = df.append(tmpdf)

        gp = df.groupby(['o'])
        means = gp.mean()
        error = gp.std()
        plt.grid(b=True, which='both', axis='y', linestyle=':')
        plt.errorbar(means.index, means['t'], yerr=error['t'], marker='.', fmt='o', mfc='b', capthick=2, ms=7, label=percentage)


    plt.savefig('C:/Users/Andres/ss-final/graphs/time.png')
    #handles, labels = get_handles_and_labels_for_sorted_legend()
    #legend = plt.legend(handles, labels, loc='best', title='% de partículas')
    #plt.setp(legend.get_title(), multialignment='center')
    plt.show()
    plt.close()

def plot_survivors_with_mean_and_error(orientation_percentages, N=20):
    plt.figure()
    plt.ylabel('Porcentaje de particulas vivas [%]')
    plt.xlabel('Distribucion de orientacion inicial en sentido horario [%]')

    df = pd.DataFrame(orientation_percentages)

    for percentage in orientation_percentages:
        data_frames = [pd.read_csv(data_file) for data_file in get_data_files(str(percentage))]
        survivors = [df.iloc[-1].n*100.0/20 for df in data_frames]
        tmpdf = pd.DataFrame(survivors, columns=['s_rate'])
        tmpdf['o'] = percentage

        df = df.append(tmpdf)

        gp = df.groupby(['o'])
        means = gp.mean()
        error = gp.std()
        plt.grid(b=True, which='both', axis='y', linestyle=':')
        plt.errorbar(means.index, means['s_rate'], yerr=error['s_rate'], marker='.', fmt='o', mfc='b', capthick=2, ms=7, label=percentage)


    #plt.savefig('C:/Users/Andres/ss-final/graphs/time.png')
    #handles, labels = get_handles_and_labels_for_sorted_legend()
    #legend = plt.legend(handles, labels, loc='best', title='% de partículas')
    #plt.setp(legend.get_title(), multialignment='center')
    plt.show()
    plt.close()

def plot_orientation_in_time_with_mean_and_error(orientation_percentages, N=20):
    plt.figure()
    plt.ylabel('Porcentaje de particulas orbitando en sentido horario [%]')
    plt.xlabel('Tiempo [s]')

    for percentage in orientation_percentages:
        data_frames = [pd.read_csv(data_file) for data_file in get_data_files(str(percentage))]

        df = pd.concat(data_frames)
        df = df[df.index % 15 == 0]
        df['clock_rate'] = 100*df['clock']/df['n']
        print(df)
        gp = df.groupby(['t'])
        means = gp.mean()
        error = gp.std()
        print(means)
        print(error)
        plt.grid(b=True, which='both', axis='y', linestyle=':')
        plt.errorbar(means.index, means['clock_rate'], yerr=error['clock_rate'], marker='.', fmt='o', mfc='b', capthick=2, ms=7, label=percentage)


    #plt.savefig('C:/Users/Andres/ss-final/graphs/time.png')
    handles, labels = get_handles_and_labels_for_sorted_legend()
    legend = plt.legend(handles, labels, loc='best', title='% de partículas')
    plt.setp(legend.get_title(), multialignment='center')
    plt.show()
    plt.close()

# plot_results('Energía total [J]', plot_total_energy)
# plot_results('Energía cinética [J]', plot_kinetic_energy)
# plot_results('% de partículas en sentido horario', plot_clockwise_particles)
# plot_results('Colisiones', plot_collisions)
# plot_results('Cantidad de partículas', plot_count)
#plot_time_with_mean_and_error(orientation_percentages=[50, 70, 90])
plot_survivors_with_mean_and_error([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
#plot_orientation_in_time_with_mean_and_error([30, 60, 80])
# plot_survivors(N=50)
