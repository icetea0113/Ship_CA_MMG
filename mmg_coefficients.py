#!/usr/bin/env python
# -*- coding: utf-8 -*-
import dataclasses
from typing import List

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

import scipy.interpolate
from scipy.misc import derivative

ρ = 1025.0

@dataclasses.dataclass
class ManeuveringParams:
    X_u_dot:float
    Y_v_dot:float
    N_r_dot:float
    X_vv:float
    X_rr:float
    X_vr:float
    Y_v:float
    Y_r:float
    Y_vvv:float
    Y_vvr:float
    Y_vrr:float
    Y_rrr:float
    N_v:float
    N_r:float
    N_vvv:float
    N_vvr:float
    N_vrr:float
    N_rrr:float
    R_0:float = -1
    X_uu:float = -1
    X_vvvv:float = 0

@dataclasses.dataclass
class PrincipalDimensions:
    scale:float
    Lpp:float
    B:float
    d:float
    displacement:float
    x_G:float
    C_b:float
    D_p:float
    H_r:float
    A_r:float
    t_p:float
    k_0:float
    k_1:float
    k_2:float
    κ:float
    ϵ:float
    x_P:float
    x_r:float
    x_h:float
    a_H:float
    f_a:float
    t_R:float
    w_P0:float
    l_R:float
    γ_R_minus:float
    γ_R_plus:float
    C_1:float = -1
    C_2_minus:float = -1
    C_2_plus:float = -1
    η:float = -1
    I_zz:float = -1

@dataclasses.dataclass
class Ship:
    def __init__(self, name, principal_dimensions:PrincipalDimensions, maneuvering_params:ManeuveringParams):
        self.principal_dimensions = principal_dimensions
        self.maneuvering_params = maneuvering_params
        self.name = name
        time: List[float] = dataclasses.field(default_factory=list)
        u: List[float] = dataclasses.field(default_factory=list)
        v: List[float] = dataclasses.field(default_factory=list)
        r: List[float] = dataclasses.field(default_factory=list)
        x: List[float] = dataclasses.field(default_factory=list)
        y: List[float] = dataclasses.field(default_factory=list)
        psi: List[float] = dataclasses.field(default_factory=list)
        δ: List[float] = dataclasses.field(default_factory=list)
        npm: List[float] = dataclasses.field(default_factory=list)

        self.principal_dimensions.I_zz = ρ * self.principal_dimensions.displacement * ((0.25 * self.principal_dimensions.Lpp) ** 2)
        self.principal_dimensions.x_h *= self.principal_dimensions.Lpp
        self.principal_dimensions.x_r *= self.principal_dimensions.Lpp
        self.F_N = []

    def register_simulation_result(self,time: List[float],u: List[float],v: List[float],r: List[float],x: List[float],y: List[float],psi: List[float]):
        self.time = time
        self.u = u
        self.v = v
        self.r = r
        self.x = x
        self.y = y
        self.psi = psi
    
    def load_simulation_result(self,time: List[float],u: List[float],v: List[float],r: List[float],x0: float = 0.0,y0: float = 0.0,psi0: float = 0.0,
    ):
        x = [x0]
        y = [y0]
        psi = [psi0]
        for i, (ut, vt, rt) in enumerate(zip(u, v, r)):
            if i > 0:
                dt = time[i] - time[i - 1]
                x.append(x[-1] + (ut * np.cos(psi[-1]) - vt * np.sin(psi[-1])) * dt)
                y.append(y[-1] + (ut * np.sin(psi[-1]) + vt * np.cos(psi[-1])) * dt)
                psi.append(psi[-1] + rt * dt)

        # Register
        self.time = time
        self.u = u
        self.v = v
        self.r = r
        self.x = x
        self.y = y
        self.psi = psi
    
    def draw_xy_trajectory(
        self,
        dimensionless: bool = False,
        aspect_equal: bool = True,
        num: int or str = None,
        figsize: List[float] = [6.4, 4.8],
        dpi: float = 100.0,
        fmt: str = None,
        facecolor: str = None,
        edgecolor: str = None,
        frameon: bool = True,
        FigureClass: matplotlib.figure.Figure = matplotlib.figure.Figure,
        clear: bool = False,
        tight_layout: bool = False,
        constrained_layout: bool = False,
        save_fig_path: str = None,
        **kwargs
    ) -> plt.Figure:
        fig = plt.figure(
                num=num,
                figsize=figsize,
                dpi=dpi,
                facecolor=facecolor,
                edgecolor=edgecolor,
                frameon=frameon,
                FigureClass=FigureClass,
                clear=clear,
                # tight_layout=tight_layout,
                constrained_layout=constrained_layout,
            )
        if dimensionless:
            if fmt is None:
                plt.plot(np.array(self.x) / self.L, np.array(self.y) / self.L, **kwargs)
            else:
                plt.plot(
                    np.array(self.x) / self.L, np.array(self.y) / self.L, fmt, **kwargs
                )
            plt.xlabel(r"$x/L$")
            plt.ylabel(r"$y/L$")
        else:
            plt.plot(self.y / self.principal_dimensions.Lpp , self.x / self.principal_dimensions.Lpp, color="limegreen" )
            plt.xlabel(r"$y$")
            plt.ylabel(r"$x$")
            plt.xlim(-2,4)
            plt.ylim(-2,4)
            plt.plot()
            if aspect_equal:
                plt.gca().set_aspect("equal")
        if save_fig_path is not None:
            save_fig_path = "result/" + save_fig_path
            plt.savefig(save_fig_path, dpi=300)
        plt.close()
        return fig
    def draw_chart(
        self,
        x_index: str,
        y_index: str,
        xlabel: str = None,
        ylabel: str = None,
        num: int or str = None,
        figsize: List[float] = [6.4, 4.8],
        dpi: float = 100.0,
        facecolor: str = None,
        edgecolor: str = None,
        frameon: bool = True,
        FigureClass: matplotlib.figure.Figure = matplotlib.figure.Figure,
        clear: bool = False,
        fmt: str = None,
        save_fig_path: str = None,
        x_lim : list[float] = [],
        y_lim : list[float] = [],
        **kwargs
    ) -> plt.Figure:
        target_x = None
        if x_index == "time":
            target_x = self.time
        elif x_index == "u":
            target_x = self.u
        elif x_index == "v":
            target_x = self.v
        elif x_index == "r":
            target_x = self.r * 180 / np.pi
        elif x_index == "x":
            target_x = self.x
        elif x_index == "y":
            target_x = self.y
        elif x_index == "psi":
            target_x = self.psi
        elif x_index == "delta":
            target_x = self.δ
        elif x_index == "δ":
            target_x = self.δ
        elif x_index == "npm":
            target_x = self.npm
        if target_x is None:
            raise Exception(
                "`x_index` is not good. Please set `x_index` from ["
                "time"
                ", "
                " u"
                ", "
                " v"
                ", "
                " r"
                ", "
                " x"
                ", "
                " y"
                ", "
                " psi"
                ", "
                " delta"
                ", "
                " δ"
                ", "
                " npm"
                "]"
            )

        target_y = None
        if y_index == "time":
            target_y = self.time
        elif y_index == "u":
            target_y = self.u
        elif y_index == "v":
            target_y = self.v
        elif y_index == "r":
            target_y = [i * 180 / np.pi for i in self.r]
        elif y_index == "x":
            target_y = self.x
        elif y_index == "y":
            target_y = self.y
        elif y_index == "psi":
            target_y = self.psi
        elif y_index == "delta":
            target_y = self.δ
        elif y_index == "δ":
            target_y = self.δ
        elif y_index == "npm":
            target_y = self.npm
        elif y_index == "beta":
            target_y = np.arctan2(-1 * self.v, self.u) * 180 / np.pi
        elif y_index == "F_N":
            target_y = self.F_N
        if target_y is None:
            raise Exception(
                "`y_index` is not good. Please set `y_index` from ["
                "time"
                ", "
                " u"
                ", "
                " v"
                ", "
                " r"
                ", "
                " x"
                ", "
                " y"
                ", "
                " psi"
                ", "
                " delta"
                ", "
                " δ"
                ", "
                " npm"
                "]"
                "]"
            )
        fig = plt.figure(
            num=num,
            figsize=figsize,
            dpi=dpi,
            facecolor=facecolor,
            edgecolor=edgecolor,
            frameon=frameon,
            FigureClass=FigureClass,
            clear=clear,
        )
        if xlabel is not None:
            plt.xlabel(xlabel)
        if ylabel is not None:
            plt.ylabel(ylabel)
        if fmt is None:
            plt.plot(target_x, target_y, color="limegreen",**kwargs)
        else:
            plt.plot(target_x, target_y, fmt,color="limegreen", **kwargs)
        # if len(x_lim) != 0:
        #     plt.xlim(x_lim[0],x_lim[1])
        #     plt.ylim(y_lim[0],y_lim[1])
        if save_fig_path is not None:
            save_fig_path = "result/" + save_fig_path
            plt.savefig(save_fig_path)
        plt.close()

        return fig
    
    def __str__(self):
        return "Ship type name: {self.name}".format(self=self)

S175 = Ship(
    "S175",
    principal_dimensions=PrincipalDimensions(
        scale = 1/50, #
        Lpp = 3.5, #
        B = 0.508, #
        d = 0.19, #
        displacement = 193.57 / 1025, #
        x_G = -0.051, #
        C_b = 0.572, #
        D_p = 0.1301, #
        H_r = 0.154, #
        A_r = 0.0130, #
        t_p=0.175, #
        k_0 = 0.2932, #
        k_1 = -0.1971, #
        k_2 = -0.0481, #
        η = 0.1301/0.154, #
        κ = 0.631, #
        ϵ = 0.921, #
        x_P = -0.5, #
        x_r = -0.5,  
        x_h = -0.48, #
        a_H = 0.237, #
        f_a = 2.75, #
        t_R = 0.29, #
        w_P0 = 0.1684, #
        l_R = -1.0, #
        γ_R_minus = 0.193, #
        γ_R_plus = 0.088, #
        C_1 = -8.0, #
    ),
    maneuvering_params=ManeuveringParams(
        X_u_dot = 0.0044, #
        X_uu = 0.01563, #
        Y_v_dot = 0.1299, #
        N_r_dot = 0.0077, #
        X_vv = -0.0711, #
        X_rr = 0.0037, #
        X_vr = -0.0726 + 0.1299, #
        Y_v = 0.2137, #
        Y_r = 0.0402 + 0.0044, #
        Y_vvv = 2.008, #
        Y_vvr = 0.3942, #
        Y_vrr = 0.7461, #
        Y_rrr = 0.0326, #
        N_v = 0.0710, #
        N_r = -0.0409, #
        N_vvv = -0.0275, #
        N_vvr = -0.7811, #
        N_vrr = -0.0287, #
        N_rrr = -0.0422 #
    )
)

KVLCC2 = Ship(
    "KVLCC2",
    principal_dimensions=PrincipalDimensions(
        scale = 1/45.7,
        Lpp = 7.00,
        B = 1.27,
        d = 0.46,
        displacement = 3.27,
        x_G = 0.25,
        C_b = 0.810,
        D_p = 0.216,
        H_r = 0.345,
        A_r = 0.0539,
        t_p = 0.220,
        k_0 = 0.2931,
        k_1 = -0.2753,
        k_2 = -0.1385,
        C_1 = 2.0,
        C_2_minus = 1.1,
        C_2_plus = 1.6,
        η = 0.216/0.345,
        κ = 0.5,
        ϵ = 1.09,
        x_P = -0.48,
        x_r = -0.5,
        x_h = -0.464,
        a_H = 0.312,
        f_a = 2.747,
        t_R = 0.387,
        w_P0 = 0.4,
        l_R = -0.710,
        γ_R_minus = 0.395,
        γ_R_plus = 0.640
    ),
    maneuvering_params=ManeuveringParams(
        X_u_dot = 0.022,
        Y_v_dot = 0.223,
        N_r_dot = 0.011,
        X_vv = -0.040,
        X_rr = 0.011,
        X_vr = 0.002,
        X_vvvv = 0.771,
        Y_v = -0.315,
        Y_r = 0.083,
        Y_vvv = -1.607,
        Y_vvr = 0.379,
        Y_vrr = -0.391,
        Y_rrr = 0.008,
        N_v = -0.137,
        N_r = -0.049,
        N_vvv = -0.030,
        N_vvr = -0.294,
        N_vrr = 0.055,
        R_0 = 0.022,
        N_rrr = -0.013
    )
)