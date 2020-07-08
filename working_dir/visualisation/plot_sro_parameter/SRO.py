import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("CSV_32x32x10-xIn=0.5_T=773K_SRO_X.csv", header=None)
df.columns = ['Number of MC Steps', 'SRO 1', 'SRO 2', 'SRO 3', 'SRO 4']

# Return Line plots
cols_list = df.columns[1:]
df.plot(x=df.columns[0], y=df.columns[1:])
plt.title("SRO Parameter vs. Time in MC Steps")
plt.ylabel('SRO Parameter')
plt.show()

# Returns scatter plots - remove ax1 variables assignments and ax-arg to print 4 plots

ax1 = df.plot(kind='scatter', x='Number of MC Steps', y='SRO 1', color='red', s=10)
ax2 = df.plot(kind='scatter', x='Number of MC Steps', y='SRO 2', color='blue', ax=ax1, s=10)
ax3 = df.plot(kind='scatter', x='Number of MC Steps', y='SRO 3', color='green', ax=ax1, s=10)
ax4 = df.plot(kind='scatter', x='Number of MC Steps', y='SRO 4', color='purple', ax=ax1, s=10)
plt.title('SRO Parameter as a function of MC timesteps')
plt.ylabel('SRO Parameter')
plt.ylim(0, 1)
plt.xlim(-100000, 5000000)
plt.xlabel('time (MC Steps)')

plt.show()
