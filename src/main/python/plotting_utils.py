import os
import pandas as pd
import matplotlib.pyplot as plt


JAR_PATH = '/home/lucas/Documents/ss-final/target/ss-final-1.0-jar-with-dependencies.jar'

# Current directory
curr_dir = os.getcwd()


def get_data_folder(sub_folder=''):
    # ../../../data/ + sub_folder
    return os.path.abspath(os.path.join(curr_dir, os.pardir, os.pardir, os.pardir, 'data/' + sub_folder))


def get_data_files(sub_folder):
    data_folder = get_data_folder(sub_folder)
    return [os.path.join(data_folder, file_name) for file_name in os.listdir(data_folder)
            if file_name.endswith('.data')]


# Copy-pasted function, not sure what is going on but it sorts the legend labels in decreasing order
def get_handles_and_labels_for_sorted_legend():
    handles, labels = plt.gca().get_legend_handles_labels()
    return zip(*[(handles[i], labels[i]) for i in reversed(sorted(range(len(handles)),
                                                                  key=lambda k: list(map(int, labels))[k]))])


def plot_mean_and_error(var_params, sub_folder_prefix):
    for var_param in var_params:
        data_frames = [pd.read_csv(data_file).iloc[::30, :] for data_file in get_data_files(sub_folder_prefix +
                                                                                            str(var_param))]
        # To plot only up to t = 1500s
        data_frames = [df[df.t <= 1500] for df in data_frames]
        gp = pd.concat(data_frames).groupby(['t'])
        means = gp.mean()
        error = gp.std()
        plt.grid(b=True, which='both', axis='y', linestyle=':')
        plt.errorbar(means.index, means['n'], yerr=error['n'], marker='.', fmt='o', mfc='c', capthick=2, ms=7,
                     label=var_param)


def average(lst):
    return sum(lst) / len(lst)