
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse


# parse arguments
parser = argparse.ArgumentParser(description='Turn dominant_processes output into L-t preview. --fraction of 0 and --digits of at least 5 passed to dominant_processes are recommended.')
parser._action_groups.pop()
requiredArgs = parser.add_argument_group('required arguments')
requiredArgs.add_argument('--input', type=str, help='dominant_processes output', required=True)
optionalArgs = parser.add_argument_group('optional arguments')
optionalArgs.add_argument('--output', type=str, help='output L_nu-t diagram', default='dominant_processes.pdf')
optionalArgs.add_argument('--time_range', type=float, nargs=2, help='time[years] range to plot', default=None)
optionalArgs.add_argument('--plot_luminosity', action='store_true', help='plot luminosity instead of percentage of each process')
optionalArgs.add_argument('--luminosity_range', type=float, nargs=2, help='luminosity[erg/s] range to plot (if --plot_luminosity is set)', default=None)

args = parser.parse_args()

# read data
df = pd.read_csv(args.input, sep='\s+', skiprows=1)

# extract names and percentages of dominant processes from the 'Dominant_processes' column
def extract_names_and_percentages(s):
    strings = s.split(',')
    names = []
    percentages = []
    for string in strings:
        percentage, name = string.split(']')
        percentage = float(percentage.strip('[').strip('%'))
        # only extract processes that actually contribute to the luminosity
        if percentage == 0:
            continue 
        names.append(name)
        percentages.append(percentage)
    return names, percentages

# determine all unique process names in the dataframe, assign them a color and a unique label for the legend
unique_names = set()
for s in df['Dominant_processes']:
    names, percentages = extract_names_and_percentages(s)
    unique_names.update(names)
unique_names = sorted(unique_names)[::-1] 
colors = plt.get_cmap('viridis', len(unique_names))
name_to_color = {name: colors(i) for i, name in enumerate(unique_names)}

# plot dominant processes
if not args.plot_luminosity:
    for i in range(len(df) - 1):
        current_time, future_time = df.loc[i, 't[years]'], df.loc[i + 1, 't[years]']
        names, percentages = extract_names_and_percentages(df.loc[i, 'Dominant_processes'])
        # reorder names and percentages according to the order of unique_names to ensure consistent coloring
        names, percentages = zip(*sorted(zip(names, percentages), key=lambda x: unique_names.index(x[0])))
        cumulative_percentages = np.cumsum(percentages)
        cumulative_percentages = np.insert(cumulative_percentages, 0, 0) 
        for j in range(len(names)):
            plt.fill_between([current_time, future_time], [cumulative_percentages[j], cumulative_percentages[j]], [cumulative_percentages[j + 1], cumulative_percentages[j + 1]], color=name_to_color[names[j]])
    plt.ylabel(r'$L^{\infty}_{proc}/L^{\infty}_{\nu}$ [%]')
else:
    # plot luminosity of each process as a function of time
    luminosities = {name : [] for name in unique_names}
    for i in range(len(df)):
        names, percentages = extract_names_and_percentages(df.loc[i, 'Dominant_processes'])
        total_luminosity = df.loc[i, 'L^inf_nu[erg/s]']
        for name, percentage in zip(names, percentages):
            luminosities[name].append(percentage / 100.0 * total_luminosity)
        for name in unique_names:
            if name not in names:
                luminosities[name].append(0.0)
    for name in unique_names:
        plt.plot(df['t[years]'], luminosities[name], color=name_to_color[name])
    plt.yscale('log')
    plt.ylabel(r'$L^{\infty}_{\nu}$ [erg/s]')
    if args.luminosity_range:
        plt.ylim(args.luminosity_range)
if args.time_range:
    plt.xlim(args.time_range)
plt.xlabel('t [years]')
plt.xscale('log')
# create legend with unique labels and colors
handles = [plt.Line2D([0], [0], color=name_to_color[name], label=name) for name in unique_names]
plt.legend(handles=handles, title='Dominant processes')

plt.savefig(args.output)