
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Process m_r_diagram raw output into M-R preview.')
parser._action_groups.pop()
requiredArgs = parser.add_argument_group('required arguments')
requiredArgs.add_argument('--input', type=str, help='tov cache to process', required=True)
optionalArgs = parser.add_argument_group('optional arguments')
optionalArgs.add_argument('--output', type=str, help='output M-R diagram', default='M-R-diagram.pdf')
args = parser.parse_args()

# read data
df = pd.read_csv(args.input, sep='\s+')

# plot M vs R
df.plot(x='R[km]', y='M[Ms]', kind='line', xlim=(0, 20), ylim=(0, 3), ylabel='M [Ms]', xlabel='R [km]')
plt.legend(['M-R diagram'])

plt.savefig(args.output)