import mmg_coefficients
import mmg_equations
import wave_data_check
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

"""
KVLCC2 Turning Test Froude Number sheet
            λ/L=0.7	    λ/L=1.0	    λ/L=1.2
Head sea	0.1	        0.079	    0.094
(χ=180°)			
Beam sea	0.138	    0.138	    0.138
(χ=270°)			

S175 Turning Test Froude Number sheet
            λ/L=0.7	    λ/L=1.0	    λ/L=1.2
Head sea	0.15	    0.15	    0.15
(χ=180°)			
Beam sea	0.15	    0.15	    0.15
(χ=270°)			

"""

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
ship_type = mmg_coefficients.S175
n_const = 10.05
# steering_rate = 1.76 * 4  # [°/s]
max_δ_rad = -35 * np.pi / 180.0  # [rad]
rudder_rate = 12 * np.pi / 180.0


# ship_type = mmg_coefficients.KVLCC2
# n_const = 10.05
# # steering_rate = 1.76 * 4  # [°/s]
# n_const = 11.8
# max_δ_rad = 35 * np.pi / 180.0  # [rad]
# rudder_rate = 11.8 * np.pi / 180.0
# u0 = 1.1767

duration = 500  # [s]

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

knot = 15
velocity = knot / 1.944

if ship_type.name == "S175":
    froude_number = velocity  / np.sqrt(9.81 * ship_type.principal_dimensions.Lpp / ship_type.principal_dimensions.scale)
    u0 = 0.879
elif ship_type.name == "KVLCC2":
    froude_number = velocity  / np.sqrt(9.81 * ship_type.principal_dimensions.Lpp / ship_type.principal_dimensions.scale)
    u0 = 15/1.9441*np.sqrt(7/320)

wave_frequency_KVLCC2 = [0.3 + 0.2 * i for i in range(0, 7)]
wave_frequency_S175 = [0.35 + 0.05 * i for i in range(0, 11)]

if ship_type.name == "S175":
    wave_frequency = 1.2
elif ship_type.name == "KVLCC2":
    wave_frequency = 1.2
wave_angle = 270
wave_amplitude = 0.01 * ship_type.principal_dimensions.Lpp

npm_list = np.array([n_const for i in range(sampling)])
wave_data_check.main(ship_type)
Wave_function_var = [wave_frequency, wave_angle, wave_amplitude]
sol = mmg_equations.simulate(ship_type, time_list,δ_rad_list,npm_list, state="wave", Wave_function_var=Wave_function_var, u0 = u0, v0=0.0, r0=0.0, method="RK45")


if sol == 0:
    Exception("Simulation failed")

result = sol.sol(time_list)

mmg_equations.calculate_drift_data()
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
