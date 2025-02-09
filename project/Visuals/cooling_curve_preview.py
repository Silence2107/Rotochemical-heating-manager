
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Turn cooling_curve (or identically structured) output into T(t) preview.')
parser._action_groups.pop()
requiredArgs = parser.add_argument_group('required arguments')
requiredArgs.add_argument(
    '--input', type=str, help='cooling t/T tabular to process', required=True)
optionalArgs = parser.add_argument_group('optional arguments')
optionalArgs.add_argument('--output', type=str, help='output t/T preview', default='Cooling.pdf')
optionalArgs.add_argument('--time_range', type=float, nargs=2, help='time[years] range to plot', default=None)
optionalArgs.add_argument('--temperature_range', type=float,
                          nargs=2, help='temperature[K] range to plot', default=None)
args = parser.parse_args()

# read data
df = pd.read_csv(args.input, sep='\s+', skiprows=1)

# plot T vs t
df.plot(x='t[years]', y='Te^inf[K]', kind='line', ylabel='Te[K], redshifted', xlabel='t[years]')
if args.time_range is not None:
    plt.xlim(args.time_range)
if args.temperature_range is not None:
    plt.ylim(args.temperature_range)
plt.legend(['RHM cooling curve'])
plt.yscale('log')
plt.xscale('log')
plt.savefig(args.output)