import pandas as pd
import matplotlib.pyplot as plt

data1 = pd.read_csv('Yasukawa_cal.csv')
data2 = pd.read_csv('Yasukawa_exp.csv')

plt.figure(figsize=(10, 6)) 

plt.plot(data1['x'], data1['Curve1'], label='Yasukawa_cal')

plt.plot(data2['x'], data2['Curve1'], label='Yasukawa_exp')

plt.title('S175 zigzag')
plt.xlabel('X Axis Label')
plt.ylabel('Y Axis Label')

print(min(data1['Curve1']), max(data1['Curve1']))
print(min(data2['Curve1']), max(data2['Curve1']))
plt.legend()

plt.show()
