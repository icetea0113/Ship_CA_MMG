import pandas as pd
import matplotlib.pyplot as plt
import mmg_equations
import mmg_coefficients
import numpy as np

knot = 6
velocity = knot / 1.944

ship_type = mmg_coefficients.KVLCC2

if ship_type.name == "S175":
    froude_number = velocity  / np.sqrt(9.81 * ship_type.principal_dimensions.Lpp / ship_type.principal_dimensions.scale)
elif ship_type.name == "KVLCC2":
    froude_number = velocity  / np.sqrt(9.81 * ship_type.principal_dimensions.Lpp / ship_type.principal_dimensions.scale)
v = froude_number * (np.sqrt(9.81 * ship_type.principal_dimensions.Lpp))

wave_frequency_KVLCC2 = [0.3 + 0.2 * i for i in range(0, 7)]
wave_frequency_S175 = [0.35 + 0.05 * i for i in range(0, 11)]
wave_frequency_S175 = [0.5, 0.7, 1.0, 1.2, 1.5]

if ship_type.name == "S175":
    wave_frequency = [i * ship_type.principal_dimensions.Lpp for i in wave_frequency_S175]
    entrance_run_angle = 52 * np.pi/180
elif ship_type.name == "KVLCC2":
    wave_frequency = [i * ship_type.principal_dimensions.Lpp for i in wave_frequency_KVLCC2]
    entrance_run_angle = 52 * np.pi/180
wave_angle = 30 * np.pi/180
wave_amplitude = 0.02 * ship_type.principal_dimensions.Lpp

wave_force = []
wave_force_dimensionless = []
for i in wave_frequency:
    X_W, Y_W, Z_W = mmg_equations.WaveForceEquation(ship_type,v, 0, 0, entrance_run_angle, i, wave_angle, wave_amplitude)
    wave_force.append(list([X_W, Y_W, Z_W]))

    X_W = X_W / (1025 * 9.81 * ship_type.principal_dimensions.B ** 2 * wave_amplitude ** 2 / ship_type.principal_dimensions.Lpp)
    Y_W = Y_W / (1025 * 9.81 * ship_type.principal_dimensions.B ** 2 * wave_amplitude ** 2 / ship_type.principal_dimensions.Lpp)
    Z_W = Z_W / (1025 * 9.81 * ship_type.principal_dimensions.B ** 3 * wave_amplitude ** 2 / ship_type.principal_dimensions.Lpp)
    wave_force_dimensionless.append(list([X_W, Y_W, Z_W]))

X_W_dimensionless = [i[0] for i in wave_force_dimensionless]
Y_W_dimensionless = [i[1] for i in wave_force_dimensionless]
Z_W_dimensionless = [i[2] for i in wave_force_dimensionless]
print(X_W_dimensionless)
plt.plot(wave_frequency_KVLCC2, Y_W_dimensionless)
plt.xlim(0,2)
plt.ylim(-10, 20)
plt.grid(False)
plt.show()