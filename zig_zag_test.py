import mmg_coefficients
import mmg_equations
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd

def main():
    ship_type = mmg_coefficients.Ship("empty", mmg_coefficients.PrincipalDimensions(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), mmg_coefficients.ManeuveringParams(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0))

    print("Excute Zig Zag Test")

    duration = 100
    num_of_sampling = 20000
    time_list = np.linspace(0.00, duration, num_of_sampling)

    print(
"**************************\n\
1. S175 Container Ship\n\
2. KVLCC2 Tanker\n\
**************************")
    
    ship_type_str = input("Enter ship type number: ")
    if ship_type_str == "1":
        ship_type = mmg_coefficients.S175
        target_δ_rad = 10.0 * np.pi / 180.0
        target_ψ_rad_deviation = 10 * np.pi / 180.0
        n_const = 7.54 # [rpm]
        npm_list = np.array([n_const for i in range(num_of_sampling)])
        u0 = 0.879
        δ_rad_rate = 12.0* np.pi / 180

    elif ship_type_str == "2":
        ship_type = mmg_coefficients.KVLCC2
        target_δ_rad = 10.0 * np.pi / 180.0
        target_ψ_rad_deviation = 10 * np.pi / 180.0
        n_const = 11.8  # [rpm]
        npm_list = np.array([n_const for i in range(num_of_sampling)])
        u0 = 1.1767
        δ_rad_rate = 11.8* np.pi / 180

    # print("************ debuging mode ************")
    # ship_type = mmg_coefficients.KVLCC2

    # target_δ_rad = 20.0 * np.pi / 180.0
    # target_ψ_rad_deviation = 20 * np.pi / 180.0
    # duration = 80
    # num_of_sampling = 8000
    # time_list = np.linspace(0.00, duration, num_of_sampling)
    # n_const = 11.8  # [rpm]
    # npm_list = np.array([n_const for i in range(num_of_sampling)])
    # u0 = 1.1767


    δ_list, u_list, v_list, r_list, x_list, y_list, ψ_list = mmg_equations.zigzag_test_mmg_3dof(
        ship_type,
        target_δ_rad,
        target_ψ_rad_deviation,
        time_list,
        npm_list,
        u0=u0,
        δ_rad_rate=δ_rad_rate,
    )

    if ship_type_str == "2":
        ship = mmg_coefficients.KVLCC2
        ship.load_simulation_result(time_list, u_list, v_list, r_list)
        ship.δ = δ_list
        ship.npm = npm_list

        ship2 = mmg_coefficients.KVLCC2
        ship2.register_simulation_result(time_list, u_list, v_list, r_list, x_list, y_list, ψ_list)
        ship2.δ = δ_list
        ship2.npm = npm_list

    elif ship_type_str == "1":
        ship = mmg_coefficients.S175
        ship.load_simulation_result(time_list, u_list, v_list, r_list)
        ship.δ = δ_list
        ship.npm = npm_list

        ship2 = mmg_coefficients.S175
        ship2.register_simulation_result(time_list, u_list, v_list, r_list, x_list, y_list, ψ_list)
        ship2.δ = δ_list
        ship2.npm = npm_list

    plt.plot(time_list, list(map(lambda psi: psi * 180 / np.pi, ship.psi)), color="limegreen")
    plt.plot(time_list, list(map(lambda δ: δ * 180 / np.pi, ship.δ)),color="limegreen")
    # plt.xlim(0, 200)
    # plt.ylim(-30, 30)
    print(max(ship.psi[0:40000])* 180 / np.pi, min(ship.psi[0:70000])* 180 / np.pi)

    ship2.draw_chart(
        "time",
        "u",
        xlabel="time [sec]",
        ylabel=r"$u$" + " [m/s]",
        save_fig_path="zigzag_u_time_series.png",
        x_lim=[0,100,10],
        y_lim=[0.6,1.4,0.2]
    )

    # ship2.draw_chart(
    #     "time",
    #     "F_N",
    #     xlabel="time [sec]",
    #     ylabel=r"$v$" + " [m/s]",
    #     save_fig_path="zigzag_v_time_series.png",
    #     x_lim=[0,200,10],
    #     y_lim=[-0.4,0.4,0.2]
    # )

    ship2.draw_chart(
        "time",
        "r",
        xlabel="time [sec]",
        ylabel=r"$r$" + " [m/s]",
        save_fig_path="zigzag_r_time_series.png",
        # x_lim=[0,100,10],
        # y_lim=[-4,4,2]
    )
    df = pd.DataFrame(time_list, columns=["time"])
    df["angle"] = ship.F_N
    df.to_csv("zigzag_only_uR.csv", index=False)

    plt.show()

if __name__ == "zig_zag_test":
    print("Zig Zag Test is ready")
elif __name__ == "main":
    print("main~!")