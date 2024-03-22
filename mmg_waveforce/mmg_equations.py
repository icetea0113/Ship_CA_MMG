import mmg_coefficients as mmg
from typing import List
import numpy as np

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.misc import derivative
import scipy.special as sp

ρ = 1025.0
Fn_global = []

def simulate_mmg_3dof(
    ship: mmg.Ship,
    time_list: List[float],
    δ_list: List[float],
    npm_list: List[float],
    u0: float = 0.0,
    v0: float = 0.0,
    r0: float = 0.0,
    x0: float = 0.0,
    y0: float = 0.0,
    ψ0: float = 0.0,
    method: str = "RK45",
    t_eval=None,
    events=None,
    vectorized=False,
    **options
):
    return simulate(
        ship,
        time_list=time_list,
        δ_list=δ_list,
        npm_list=npm_list,
        u0=u0,
        v0=v0,
        r0=r0,
        x0=x0,
        y0=y0,
        ψ0=ψ0,
        method=method,
        t_eval=t_eval,
        events=events,
        vectorized=vectorized,
        **options
    )
def simulate(ship: mmg.Ship, time_list: List[float], δ_list: List[float], npm_list: List[float], Wave_function_var = [], u0: float = 0.0, v0: float = 0.0, r0: float = 0.0, x0: float = 0.0, y0: float = 0.0, ψ0: float = 0.0, method: str = "RK45", t_eval=None, events=None, vectorized=False, state = "No-wave", **options):

    spl_δ = interp1d(time_list, δ_list, "cubic", fill_value="extrapolate")
    spl_npm = interp1d(time_list, npm_list, "cubic", fill_value="extrapolate")

    def MMG_3dof_eom_solve_ivp(t, X):
        u, v, r, x, y, ψ, δ, npm = X
        d_u, d_v, d_r = MotionEquation(ship, npm, u, v, r, δ, Wave_function_var, state)
        d_x = u * np.cos(ψ) - v * np.sin(ψ)
        d_y = u * np.sin(ψ) + v * np.cos(ψ)
        d_ψ = r
        d_δ = derivative(spl_δ, t)
        d_npm = derivative(spl_npm, t)
        return [d_u, d_v, d_r, d_x, d_y, d_ψ, d_δ, d_npm]

    sol = solve_ivp(
    MMG_3dof_eom_solve_ivp,
    [time_list[0], time_list[-1]],
    [u0, v0, r0, x0, y0, ψ0, δ_list[0], npm_list[0]],
    dense_output=True,
    method=method,
    t_eval=t_eval,
    events=events,
    vectorized=vectorized,
    **options
    )
    return sol
    
def MotionEquation(ship: mmg.Ship, npm, u, v, r, δ, Wave_function_var, state = "No-wave"):
    X_H, Y_H, N_H = HullForcesEquation_no_Added_mass(ship, u, v, r)
    X_P, X_R, Y_R, N_R = Propeller_RudderForcesEquation(ship, npm, u, v, r, δ)

    if state == "Wave":
        entrance_run_angle, wave_frequency, wave_angle, wave_amplitude = Wave_function_var
        X_W, Y_W, N_W = WaveForceEquation(ship, u, v, r, entrance_run_angle, wave_frequency, wave_angle, wave_amplitude)
    elif state == "No-wave":
        X_W, Y_W, N_W = 0, 0, 0
    else:
        raise ValueError("Invalid state")

    v_m = v - ship.principal_dimensions.x_G * r
    m = ship.principal_dimensions.displacement*ρ
    m_x = ship.maneuvering_params.X_u_dot * (0.5 * ρ * (ship.principal_dimensions.Lpp ** 2) * ship.principal_dimensions.d)
    m_y = ship.maneuvering_params.Y_v_dot * (0.5 * ρ * (ship.principal_dimensions.Lpp ** 2) * ship.principal_dimensions.d)
    J_x = ship.maneuvering_params.N_r_dot * (0.5 * ρ * (ship.principal_dimensions.Lpp ** 4) * ship.principal_dimensions.d)

    du = ((m + m_y)*v_m*r + ship.principal_dimensions.x_G * m * r ** 2 + X_R + X_P + X_H + X_W) / (m + m_x)

    dv_m = ((ship.principal_dimensions.x_G**2) * (m**2) * u * r
            - (N_H + N_R + N_W) * ship.principal_dimensions.x_G * m
            + ((Y_H + Y_R + Y_W) - (m + m_x) * u * r) * (ship.principal_dimensions.I_zz + J_x + (ship.principal_dimensions.x_G**2) * m)
        ) / ((ship.principal_dimensions.I_zz + J_x + (ship.principal_dimensions.x_G**2) * m) * (m + m_y) - (ship.principal_dimensions.x_G**2) * (m**2))

    dr = (N_R + N_H + N_W - ship.principal_dimensions.x_G * m * (dv_m + u * r)) / (ship.principal_dimensions.I_zz + J_x + ship.principal_dimensions.x_G ** 2 * m)

    dv = dv_m + ship.principal_dimensions.x_G * dr
    print("X_H: ", X_H, ", Y_H: ", Y_H, ", N_H: ", N_H, ", X_P: ", X_P, ", X_R: ", X_R, ", Y_R: ", Y_R, ", N_R: ", N_R, ", X_W: ", X_W, ", Y_W: ", Y_W, ", N_W: ", N_W, "\n")
    return du, dv, dr

def HullForcesEquation_no_Added_mass(ship: mmg.Ship, u, v, r):
    global ρ
    v_m = v - ship.principal_dimensions.x_G * r
    U = np.sqrt(u**2 + v_m**2)
    beta = 0.0 if U == 0.0 else np.arcsin(-v_m/U)
    v_dash = 0.0 if U == 0.0 else v_m / U
    beta_or_v_dash = v_dash if ship.name != "S175" else beta
    r_dash = 0.0 if U == 0.0 else r * ship.principal_dimensions.Lpp / U

    R = 0.5 * ρ * ship.principal_dimensions.Lpp * ship.principal_dimensions.d * (u ** 2) * ship.maneuvering_params.X_uu if ship.maneuvering_params.X_uu != -1 \
    else 0.5 * ρ * ship.principal_dimensions.Lpp * ship.principal_dimensions.d * (U ** 2) * ship.maneuvering_params.R_0

    X_H = 0.5 * ρ * ship.principal_dimensions.Lpp * ship.principal_dimensions.d * (U ** 2) * (
        ship.maneuvering_params.X_vv * (beta_or_v_dash ** 2) +
        ship.maneuvering_params.X_rr * (r_dash ** 2) +
        ship.maneuvering_params.X_vr * beta_or_v_dash * r_dash
        + ship.maneuvering_params.X_vvvv * (beta_or_v_dash ** 4)
    ) - R

    Y_H = 0.5 * ρ * ship.principal_dimensions.Lpp * ship.principal_dimensions.d * (U ** 2) * (
        ship.maneuvering_params.Y_v * beta_or_v_dash +
        ship.maneuvering_params.Y_r * r_dash +
        ship.maneuvering_params.Y_vvv * (beta_or_v_dash ** 3) +
        ship.maneuvering_params.Y_vvr * (beta_or_v_dash ** 2) * r_dash +
        ship.maneuvering_params.Y_vrr * beta_or_v_dash * (r_dash ** 2) +
        ship.maneuvering_params.Y_rrr * (r_dash ** 3)
    )

    N_H = 0.5 * ρ * ship.principal_dimensions.Lpp ** 2 * ship.principal_dimensions.d * (U ** 2) * (
        ship.maneuvering_params.N_v * beta_or_v_dash +
        ship.maneuvering_params.N_r * r_dash +
        ship.maneuvering_params.N_vvv * (beta_or_v_dash ** 3) +
        ship.maneuvering_params.N_vvr * (beta_or_v_dash ** 2) * r_dash +
        ship.maneuvering_params.N_vrr * beta_or_v_dash * (r_dash ** 2) +
        ship.maneuvering_params.N_rrr * (r_dash ** 3)
    )
    return X_H, Y_H, N_H

def Propeller_RudderForcesEquation(ship: mmg.Ship, npm, u, v, r, δ, state = 0):
    global ρ
    v_m = v - ship.principal_dimensions.x_G * r  
    U = np.sqrt(u**2 + v_m**2)
    r_dash = 0.0 if U == 0.0 else r * ship.principal_dimensions.Lpp / U

    β = 0.0 if U == 0.0 else np.arctan2(-v_m, u)
    β_p = β - ship.principal_dimensions.x_P * r_dash
    if ship.name != "S175":
        C_2 = ship.principal_dimensions.C_2_minus if β_p < 0.0 else ship.principal_dimensions.C_2_plus

    w_P = ship.principal_dimensions.w_P0 * np.exp(ship.principal_dimensions.C_1 * (β_p) ** 2) if ship.name == "S175" else \
    -1 * ((1 + (1 - np.exp(-1 * ship.principal_dimensions.C_1*np.abs(β_p))) * (C_2 - 1)) * (1 - ship.principal_dimensions.w_P0) - 1)

    J_P = 0.0 if npm == 0.0 else u*(1-w_P)/(npm * ship.principal_dimensions.D_p)
    K_T = ship.principal_dimensions.k_0 + ship.principal_dimensions.k_1 * J_P + ship.principal_dimensions.k_2 * (J_P ** 2)
    X_P = (1 - ship.principal_dimensions.t_p) * ρ * (npm ** 2) * (ship.principal_dimensions.D_p ** 4) * K_T

    β_R = β - ship.principal_dimensions.l_R * r_dash


    if ship.name != "S175":
        γ_R = ship.principal_dimensions.γ_R_minus if  β_R < 0.0 else ship.principal_dimensions.γ_R_plus
        u_R = (
                np.sqrt(ship.principal_dimensions.η * (ship.principal_dimensions.κ * ship.principal_dimensions.ϵ * 8.0 * ship.principal_dimensions.k_0 * npm**2 * ship.principal_dimensions.D_p**4 / np.pi) ** 2)
                if J_P == 0.0
                else
                u 
                * (1 - w_P)
                * ship.principal_dimensions.ϵ
                * np.sqrt(
                    ship.principal_dimensions.η * (1.0 + ship.principal_dimensions.κ * (np.sqrt(1.0 + 8.0 * K_T / (np.pi * J_P**2)) - 1)) ** 2
                    + (1 - ship.principal_dimensions.η)
                )
            )
    elif ship.name == "S175":
        γ_R = ship.principal_dimensions.γ_R_minus if  v >= 0.0 else ship.principal_dimensions.γ_R_plus

        u_R = (u 
            * (1 - w_P)
            * ship.principal_dimensions.ϵ
            * np.sqrt(1.0 + ship.principal_dimensions.κ *(8.0 * K_T / (np.pi * J_P**2)))
        )
    '''
    # - Trajectory correct but zigzag test incorrect
        # u_R = (
        #     np.sqrt(ship.principal_dimensions.η * (ship.principal_dimensions.κ * ship.principal_dimensions.ϵ * 8.0 * ship.principal_dimensions.k_0 * npm**2 * ship.principal_dimensions.D_p**4 / np.pi) ** 2)
        #     if J_P == 0.0
        #     else
        #     u 
        #     * (1 - w_P)
        #     * ship.principal_dimensions.ϵ
        #     * np.sqrt(
        #         ship.principal_dimensions.η * (1.0 + ship.principal_dimensions.κ * (np.sqrt(1.0 + 8.0 * K_T / (np.pi * J_P**2)) - 1)) ** 2
        #     ))
    '''
    v_R = U * γ_R * β_R

    a_R = δ - np.arctan2(v_R, u_R)
    U_R = np.sqrt(u_R**2 + v_R**2)

    F_N = 0.5 * ρ * ship.principal_dimensions.A_r * (U_R ** 2) * ship.principal_dimensions.f_a * np.sin(a_R)
    if state == 1:
        global Fn_global
        Fn_global.append(F_N/9.81)

    X_R = -(1-ship.principal_dimensions.t_R)*F_N*np.sin(δ)

    Y_R = -(1+ship.principal_dimensions.a_H)*F_N*np.cos(δ)

    N_R = -(ship.principal_dimensions.x_r + ship.principal_dimensions.a_H * ship.principal_dimensions.x_h) * F_N * np.cos(δ)

    return X_P, X_R, Y_R, N_R

def WaveForceEquation(ship: mmg.Ship, u, v, r, entrance_run_angle, wave_frequency, wave_angle, wave_amplitude):
    if ship.name == "KVLCC2":
        k = 2 * np.pi / wave_frequency
    elif ship.name == "S175":
        k = wave_frequency ** 2 / 9.81
        k = 2 * np.pi / wave_frequency
    d = ship.principal_dimensions.d
    b = ship.principal_dimensions.B / ship.principal_dimensions.Lpp
    I1 = sp.iv(1, 1.5*k*d)
    K1 = sp.kv(1, 1.5*k*d)
    alpha_1 = (np.pi ** 2 * I1 ** 2)/(np.pi**2 * I1**2 + K1**2)
    mu = b / np.tan(entrance_run_angle)
    # print("alpha_1: ", alpha_1, ", mu: ", mu, ", k: ", k, ", I1: ", I1, ", K1: ", K1, ", entrance_run angle: ", entrance_run_angle, ", wave_angle: ", wave_angle, "\n")

    X_W_dash = -0.5 * alpha_1 * b * (
    np.sin(entrance_run_angle + wave_angle) * abs(np.sin(entrance_run_angle + wave_angle)) +
    np.sin(entrance_run_angle - wave_angle) * abs(np.sin(entrance_run_angle - wave_angle)))
    Y_W_dash = 1 * alpha_1 * (
        (1 - mu) *  np.sin(wave_angle) * abs(np.sin(wave_angle)) +
        0.5 * mu * (np.sin(entrance_run_angle + wave_angle) * abs(np.sin(entrance_run_angle + wave_angle)) -
        np.sin(entrance_run_angle - wave_angle) * abs(np.sin(entrance_run_angle - wave_angle)))
    )
    N_W_dash = alpha_1 * b * (0.5 * mu - np.cos(entrance_run_angle) ** 2) * np.sin(wave_angle) * np.cos(wave_angle)

    X_W = 0.5 * ρ * 9.81 * ship.principal_dimensions.Lpp * wave_amplitude ** 2 * X_W_dash
    Y_W = 0.5 * ρ * 9.81 * ship.principal_dimensions.Lpp * wave_amplitude ** 2 * Y_W_dash
    N_W = 0.5 * ρ * 9.81 * ship.principal_dimensions.Lpp ** 2 * wave_amplitude ** 2 * N_W_dash

    return X_W, Y_W, N_W

def zigzag_test_mmg_3dof(
    ship: mmg.Ship,
    target_δ_rad: float,
    target_ψ_rad_deviation: float,
    time_list: List[float],
    npm_list: List[float],
    δ0: float = 0.0,
    δ_rad_rate: float = 1.0 * np.pi / 180,
    u0: float = 0.0,
    v0: float = 0.0,
    r0: float = 0.0,
    x0: float = 0.0,
    y0: float = 0.0,
    ψ0: float = 0.0,
    method: str = "RK45",
    t_eval=None,
    events=None,
    vectorized=False,
    **options
):
    print(ship)
    target_ψ_rad_deviation = np.abs(target_ψ_rad_deviation)

    final_δ_list = [0.0] * len(time_list)
    final_u_list = [0.0] * len(time_list)
    final_v_list = [0.0] * len(time_list)
    final_r_list = [0.0] * len(time_list)
    final_x_list = [0.0] * len(time_list)
    final_y_list = [0.0] * len(time_list)
    final_ψ_list = [0.0] * len(time_list)

    next_stage_index = 0
    target_δ_rad = -target_δ_rad  # for changing in while loop
    ψ = ψ0

    while next_stage_index < len(time_list):
        target_δ_rad = -target_δ_rad
        start_index = next_stage_index

        # Make delta list
        δ_list = [0.0] * (len(time_list) - start_index)
        if start_index == 0:
            δ_list[0] = δ0
            u0 = u0
            v0 = v0
            r0 = r0
            x0 = x0
            y0 = y0
        else:
            δ_list[0] = final_δ_list[start_index - 1]
            u0 = final_u_list[start_index - 1]
            v0 = final_v_list[start_index - 1]
            r0 = final_r_list[start_index - 1]
            x0 = final_x_list[start_index - 1]
            y0 = final_y_list[start_index - 1]

        for i in range(start_index + 1, len(time_list)):
            Δt = time_list[i] - time_list[i - 1]
            if target_δ_rad > 0:
                δ = δ_list[i - 1 - start_index] + δ_rad_rate * Δt
                if δ >= target_δ_rad:
                    δ = target_δ_rad
                δ_list[i - start_index] = δ
            elif target_δ_rad <= 0:
                δ = δ_list[i - 1 - start_index] - δ_rad_rate * Δt
                if δ <= target_δ_rad:
                    δ = target_δ_rad
                δ_list[i - start_index] = δ

        sol = simulate_mmg_3dof(
            ship,
            time_list[start_index:],
            δ_list,
            npm_list[start_index:],
            u0=u0,
            v0=v0,
            r0=r0,
            x0=x0,
            y0=y0,
            ψ0=ψ,
            # TODO
        )
        sim_result = sol.sol(time_list[start_index:])
        u_list = sim_result[0]
        v_list = sim_result[1]
        r_list = sim_result[2]
        x_list = sim_result[3]
        y_list = sim_result[4]
        ψ_list = sim_result[5]
        # ship = ShipObj3dof(L=basic_params.L_pp, B=basic_params.B)
        # ship.load_simulation_result(time_list, u_list, v_list, r_list, psi0=ψ)

        target_ψ_rad = ψ0 + target_ψ_rad_deviation
        if target_δ_rad < 0:
            target_ψ_rad = ψ0 - target_ψ_rad_deviation
        # ψ_list = ship.psi
        bool_ψ_list = [True if ψ < target_ψ_rad else False for ψ in ψ_list]
        if target_δ_rad < 0:
            bool_ψ_list = [True if ψ > target_ψ_rad else False for ψ in ψ_list]
        over_index_list = [i for i, flag in enumerate(bool_ψ_list) if flag is False]
        next_stage_index = len(time_list)
        if len(over_index_list) > 0:
            ψ = ψ_list[over_index_list[0]]
            next_stage_index = over_index_list[0] + start_index
            final_δ_list[start_index:next_stage_index] = δ_list[: over_index_list[0]]
            final_u_list[start_index:next_stage_index] = u_list[: over_index_list[0]]
            final_v_list[start_index:next_stage_index] = v_list[: over_index_list[0]]
            final_r_list[start_index:next_stage_index] = r_list[: over_index_list[0]]
            final_x_list[start_index:next_stage_index] = x_list[: over_index_list[0]]
            final_y_list[start_index:next_stage_index] = y_list[: over_index_list[0]]
            final_ψ_list[start_index:next_stage_index] = ψ_list[: over_index_list[0]]
        else:
            final_δ_list[start_index:next_stage_index] = δ_list
            final_u_list[start_index:next_stage_index] = u_list
            final_v_list[start_index:next_stage_index] = v_list
            final_r_list[start_index:next_stage_index] = r_list
            final_x_list[start_index:next_stage_index] = x_list
            final_y_list[start_index:next_stage_index] = y_list
            final_ψ_list[start_index:next_stage_index] = ψ_list
    for t, u, v, r, δ, npm in zip(time_list, final_u_list, final_v_list, final_r_list, final_δ_list, npm_list):
            # u, v, r, δ, npm 값에 기반해 해당 시점에서의 F_N 값을 계산
            Propeller_RudderForcesEquation(ship, npm, u, v, r, δ, 1)  # 이 함수를 수정하여 F_N을 반환하도록 함
    print(len(Fn_global))
    ship.F_N = Fn_global
    return (
        final_δ_list,
        final_u_list,
        final_v_list,
        final_r_list,
        final_x_list,
        final_y_list,
        final_ψ_list,
    )
