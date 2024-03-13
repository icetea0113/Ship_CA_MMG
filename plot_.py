import pandas as pd
import matplotlib.pyplot as plt

data1 = pd.read_csv('zigzag.csv')
data2 = pd.read_csv('zigzag_only_uR.csv')

plt.figure(figsize=(10, 6)) 

plt.plot(data1['time'], data1['angle'], label='F_N_KVVLCC2+S175')

plt.plot(data2['time'], data2['angle'], label='F_N_only_uR')

plt.title('S175 F_N')
plt.xlabel('X Axis Label')
plt.ylabel('Y Axis Label')

plt.legend()

plt.show()
