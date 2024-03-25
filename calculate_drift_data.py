from scipy.spatial import ConvexHull
import numpy as np
import matplotlib.pyplot as plt
   
angle_data = []
time_data = []
x_data = []
y_data = []
with open("ship_gps_time_angle_data.txt", 'r') as f:
    for line in f:
        data = line.split(',')
        if abs(float(data[3]) * 180 / np.pi) > 270 : 
            time_data.append(float(data[0]))
            angle_data.append(float(data[3]))
            x_data.append(float(data[2]))
            y_data.append(float(data[1]))


# # 2D 점 집합 생성
# points = np.column_stack((x_data, y_data))

# # Convex hull 찾기
# hull = ConvexHull(points)

# # 원래의 점 집합 그리기
# plt.scatter(x_data, y_data)

# # Convex hull 그리기
# for i in range(len(hull.vertices)):
#     start_vertex = hull.points[hull.vertices[i - 1]]
#     end_vertex = hull.points[hull.vertices[(i + 1) % len(hull.vertices)]]
#     plt.plot([start_vertex[0], end_vertex[0]], [start_vertex[1], end_vertex[1]], 'r-')

# # 각 점에서의 접선 그리기
# for i in range(len(hull.vertices)):
#     # 이전 점과 다음 점 찾기
#     prev_vertex = hull.points[hull.vertices[i - 1]]
#     next_vertex = hull.points[hull.vertices[(i + 1) % len(hull.vertices)]]

#     # 기울기 계산
#     slope = (next_vertex[1] - prev_vertex[1]) / (next_vertex[0] - prev_vertex[0])
#     print(f"Slope at point {hull.points[hull.vertices[i]]} is {slope}")
#     # 접선 그리기
#     x_values = np.linspace(prev_vertex[0], next_vertex[0], 100)
#     y_values = slope * (x_values - prev_vertex[0]) + prev_vertex[1]
#     plt.plot(x_values, y_values, 'g-')

# plt.show()

import matplotlib.pyplot as plt
import numpy as np

# 원래의 점 집합 그리기
plt.scatter(x_data, y_data, s=1)

# 사용자가 클릭한 두 점 얻기
print("Please click")
x = plt.ginput(4)
print("clicked", x)

# 기울기 계산
slope = (x[1][1] - x[0][1]) / (x[1][0] - x[0][0])
print(f"Slope between points {x[0]} and {x[1]} is {slope}")

print(np.sqrt((x[2][0] - x[3][0])**2 + (x[2][1] - x[3][1])**2))

# 기울기를 각도로 변환
angle_rad = np.arctan(slope)
angle_deg = np.degrees(angle_rad)

print(f"Angle in radians: {angle_rad}")
print(f"Angle in degrees: {angle_deg}")

# 접선 그리기
x_values = np.linspace(x[0][0], x[1][0], 100)
y_values = slope * (x_values - x[0][0]) + x[0][1]
plt.plot(x_values, y_values, 'g-')

plt.show()