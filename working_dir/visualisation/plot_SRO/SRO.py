import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("CSV_32x32x10-xIn=0.5_T=773K_SRO_X.csv", header=None)
df.columns = ['Number of MC Steps', 'SRO 1', 'SRO 2', 'SRO 3', 'SRO 4']
print(df.head())

ax = plt.gca
df.plot(kind='scatter', x='Number of MC Steps', y='SRO 1', color='red')
df.plot(kind='scatter', x='Number of MC Steps', y='SRO 2', color='blue')
df.plot(kind='scatter', x='Number of MC Steps', y='SRO 3', color='green')
df.plot(kind='scatter', x='Number of MC Steps', y='SRO 4', color='purple')


plt.show()
