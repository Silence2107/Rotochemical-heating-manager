
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Turn nonequilibrium_profiles output into T(r) timeline preview.')
parser._action_groups.pop()
requiredArgs = parser.add_argument_group('required arguments')
requiredArgs.add_argument(
    '--input', type=str, help='cooling r/T/t profile to process', required=True)
optionalArgs = parser.add_argument_group('optional arguments')
optionalArgs.add_argument(
    '--output', type=str, help='output profiles preview', default='CoolingProfiles.pdf')
optionalArgs.add_argument('--temperature_range', type=float,
                          nargs=2, help='temperature[K] range to plot', default=[1e6, 1e10])
optionalArgs.add_argument('--time_range', type=float,
                            nargs=2, help='time[years] range to plot', default=None)
args = parser.parse_args()

# read data, find the first row with columns
with open(args.input) as f:
    for i, line in enumerate(f):
        if not line.startswith('r[km]'):
            continue
        else:
            skiprows = i
            break

# read data
df = pd.read_csv(args.input, sep='\s+', skiprows=skiprows)

if args.time_range:
    def time_filter(column):
        return args.time_range[0] <= float(column.strip('[years]')) <= args.time_range[1]
    appropriate_columns = df.columns[1:][df.columns[1:].map(time_filter)]
else:
    appropriate_columns = df.columns[1:]

# plot T(r) for every column of interest
for column in appropriate_columns:
    plt.plot(df['r[km]'], df[column])
    # label along the line
    plt.text(df['r[km]'].iloc[0], df[column].iloc[0] * 1.05,
             column, fontdict={'size': 10, 'family': 'serif'})
plt.ylim(args.temperature_range)
plt.xlabel('r[km]')
plt.ylabel('T[K], redshifted')
plt.yscale('log')
plt.savefig(args.output)