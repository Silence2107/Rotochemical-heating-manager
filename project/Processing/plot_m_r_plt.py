
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Process m_r_diagram data into M-R preview.')
parser.add_argument('--file', type=str, help='text file to process')
parser.add_argument('--output', type=str, help='pdf output file')
args = parser.parse_args()

# read data
df = pd.read_csv(args.file, sep='\s+', skiprows=1, names=["rho [df. units]", "M [Ms]", "R [km]",])

# plot M vs R
df.plot(x='R [km]', y='M [Ms]', kind='line', title='M-R diagram', xlim=(0, 20), ylim=(0, 3), ylabel='M [Ms]', xlabel='R [km]')
plt.savefig(args.output)