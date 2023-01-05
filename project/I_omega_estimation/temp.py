
import pandas as pd
import matplotlib.pyplot as plt

inf = open('../../data/apr4_I_omega_i_lin_interp.txt', 'r')

# transform to pandas
df = pd.read_csv(inf, sep=' ')
# plot columns 1,2,4 against first one
df.plot(x=0, y=[1,2,4])
plt.ylabel('I_omega_i, s^2')
plt.show()
