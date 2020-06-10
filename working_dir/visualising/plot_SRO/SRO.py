import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("CSV_32x32x5-xIn=0.5_T=773K_SRO_1.csv", header=None)
data = data.T
print(data.head())

row = data.iloc[0]
row.plot(kind='line', legend='SRO1')

row = data.iloc[1]
row.plot(kind='line', legend='SRO2')

row = data.iloc[2]
row.plot(kind='line', legend='SRO3')

row = data.iloc[3]
row.plot(kind='line', legend='SRO4')



plt.show()
