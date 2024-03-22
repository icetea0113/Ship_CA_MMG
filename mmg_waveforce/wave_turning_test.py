import mmg_coefficients
import mmg_equations
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


ship_type = mmg_coefficients.Ship("empty", mmg_coefficients.PrincipalDimensions(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), mmg_coefficients.ManeuveringParams(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0))

print("Excute Wave Turning Test")

# print(
# "**************************\n\
# 1. S175 Container Ship\n\
# 2. KVLCC2 Tanker\n\
# **************************")

# ship_type_str = input("Enter ship type: ")
# if ship_type_str == "1":
#     ship_type = mmg_coefficients.S175
#     print(ship_type)
#     n_const = 10.05
#     max_δ_rad = -35 * np.pi / 180.0  # [rad]
#     rudder_rate = 12 * np.pi / 180.0
#     u0 = 0.879

# elif ship_type_str == "2":
#     ship_type = mmg_coefficients.KVLCC2
#     print(ship_type)
#     n_const = 11.8
#     max_δ_rad = -35 * np.pi / 180.0  # [rad]
#     rudder_rate = 11.8 * np.pi / 180.0
#     u0 = 1.1767
# else:
#     Exception("Invalid ship type")

print("************ debuging mode ************")
ship_type = mmg_coefficients.KVLCC2
n_const = 11.8
# steering_rate = 1.76 * 4  # [°/s]
max_δ_rad = 35 * np.pi / 180.0  # [rad]
rudder_rate = 11.8 * np.pi / 180.0
u0 = 1.1767

duration = 2000  # [s]

growing_rate = 10
sampling = duration * growing_rate
time_list = np.linspace(0.00, duration, sampling)
δ_rad_list = [0] * sampling
for i in range(1, len(time_list)):
    Δt = time_list[i] - time_list[i - 1]
    if max_δ_rad > 0:
        δ = δ_rad_list[i - 1] + rudder_rate * Δt
        if δ >= max_δ_rad:
            δ = max_δ_rad
        δ_rad_list[i] = δ
    elif max_δ_rad <= 0:
        δ = δ_rad_list[i - 1] - rudder_rate * Δt
        if δ <= max_δ_rad:
            δ = max_δ_rad
        δ_rad_list[i] = δ
sol = 0

knot = 6
velocity = knot / 1.944

if ship_type.name == "S175":
    froude_number = velocity  / np.sqrt(9.81 * ship_type.principal_dimensions.Lpp / ship_type.principal_dimensions.scale)
elif ship_type.name == "KVLCC2":
    froude_number = velocity  / np.sqrt(9.81 * ship_type.principal_dimensions.Lpp / ship_type.principal_dimensions.scale)
u0 = froude_number * (np.sqrt(9.81 * ship_type.principal_dimensions.Lpp))

wave_frequency_KVLCC2 = [0.3 + 0.2 * i for i in range(0, 7)]
wave_frequency_S175 = [0.35 + 0.05 * i for i in range(0, 11)]
wave_frequency_S175 = [0.5, 0.7, 1.0, 1.2, 1.5]

if ship_type.name == "S175":
    wave_frequency = [i * ship_type.principal_dimensions.Lpp for i in wave_frequency_S175]
    entrance_run_angle = 52 * np.pi/180
elif ship_type.name == "KVLCC2":
    wave_frequency = 0.7 * ship_type.principal_dimensions.Lpp
    entrance_run_angle = 52 * np.pi/180
wave_angle = 90 * np.pi/180
wave_amplitude = 0.02 * ship_type.principal_dimensions.Lpp

npm_list = np.array([n_const for i in range(sampling)])
Wave_function_var = [entrance_run_angle, wave_frequency, wave_angle, wave_amplitude]
sol = mmg_equations.simulate(ship_type, time_list,δ_rad_list,npm_list, state="No-wave", Wave_function_var=Wave_function_var, u0 = u0, v0=0.0, r0=0.0, method="RK45")


if sol == 0:
    Exception("Simulation failed")

result = sol.sol(time_list)

ship_type.register_simulation_result(time_list, result[0], result[1], result[2], result[3], result[4], result[5])

ship_type.draw_xy_trajectory(save_fig_path="trajectory.png")

ship_type.draw_chart(
    "time",
    "u",
    xlabel="time [sec]",
    ylabel=r"$u$" + " [m/s]",
    save_fig_path="u_time_series.png",
    x_lim=[0, 100, 20],
    y_lim=[0.3, 1.0, 0.1]
)

ship_type.draw_chart(
    "time",
    "v",
    xlabel="time [sec]",
    ylabel=r"$v$" + " [m/s]",
    save_fig_path="v_time_series.png",
    x_lim=[0, 100, 20],
    y_lim=[0, 0.2, 0.05]
)

ship_type.draw_chart(
    "time",
    "r",
    xlabel="time [sec]",
    ylabel=r"$r$" + " [m/s]",
    save_fig_path="r_time_series.png",
    x_lim=[0, 100, 20],
    y_lim=[-7, 0, 1]
)

ship_type.draw_chart(
    "time",
    "beta",
    xlabel="time [sec]",
    ylabel=r"$r$" + " [m/s]",
    save_fig_path="beta_time_series.png",
    x_lim=[0, 100, 20],
    y_lim=[-15, 0, 5]
)
for i in range(len(ship_type.psi)):
    if np.abs(ship_type.psi[i]) >= 90 * np.pi / 180:
        print(ship_type.x[i]/ship_type.principal_dimensions.Lpp)
        break
for i in range(len(ship_type.psi)):
    if np.abs(ship_type.psi[i]) >= 180 * np.pi / 180:
        print(ship_type.y[i]/ship_type.principal_dimensions.Lpp)
        break
    
plt.show()
