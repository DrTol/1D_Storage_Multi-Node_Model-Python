# -*- coding: utf-8 -*-
"""
Multi-Node Method - 1D Stratified Sensible Storage Unit

Developed on Sat Oct 31 12:18:27 2020
@author: Hakan İbrahim Tol, PhD

References:
[1] Kleinbach, Eberhard Markus; Performance study of one-dimensional models 
    for stratified thermal storage tank.
[2] Unrau, Cody; Numerical investigation of one-dimensional storage tank models
    and the development of analytical modelling techniques.
"""

"""
storage1d.py
Multi-Node Method - 1D Stratified Sensible Storage Unit

Method: 1D multi-node (finite volumes), with optional axial conduction and
uniform side losses. Charging mixes at inlet node; discharging pulls from
the flow direction side to preserve stratification.

Updated: Fri Aug 29 18:25:43 2025
"""

from dataclasses import dataclass
from typing import Dict
import numpy as np
import time

# --------- Small utilities ---------

def _as_scalar(x) -> float:
    """Safely convert a scalar-like value (NumPy/Python) to float."""
    if np.isscalar(x):
        return float(x)
    arr = np.asarray(x)
    if arr.ndim == 0:
        return float(arr.item())
    if arr.size == 1:
        return float(arr.reshape(()).item())
    raise ValueError("Expected a scalar-like value, got an array with size > 1.")

def _ensure_1d(arr, name: str) -> np.ndarray:
    arr = np.asarray(arr, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"'{name}' must be 1D; got shape {arr.shape}.")
    return arr

def _check_len(name: str, arr: np.ndarray, n: int):
    if len(arr) != n:
        raise ValueError(f"Length of '{name}' must be {n}, got {len(arr)}.")

# --------- Configuration ---------

@dataclass
class TankConfig:
    n_nodes: int                 # number of vertical nodes (>=2)
    height: float                # [m] total tank height
    volume: float                # [m^3] total fluid volume
    rho: float = 997.0           # [kg/m^3] water density
    cp: float = 4180.0           # [J/(kg*K)] specific heat
    ua_loss: float = 0.0         # [W/K] overall loss coefficient to ambient
    k_cond: float = 0.0          # [W/(m*K)] axial conduction (0 to disable)
    dt: float = 60.0             # [s] time step
    n_steps: int = 3600          # number of time steps

    # Port indices (0 = bottom, n_nodes-1 = top)
    charge_inlet: int = 0        # default: charge from bottom
    discharge_outlet: int = -1   # default: discharge from top (-1 normalized)

    def validate(self):
        if self.n_nodes < 2:
            raise ValueError("n_nodes must be >= 2.")
        if self.height <= 0:
            raise ValueError("height must be > 0.")
        if self.volume <= 0:
            raise ValueError("volume must be > 0.")
        if self.rho <= 0 or self.cp <= 0:
            raise ValueError("rho and cp must be > 0.")
        if self.dt <= 0:
            raise ValueError("dt must be > 0.")
        if self.n_steps < 1:
            raise ValueError("n_steps must be >= 1.")
        if not (0 <= self.charge_inlet < self.n_nodes):
            raise ValueError("charge_inlet index out of range.")
        if self.discharge_outlet == -1:
            self.discharge_outlet = self.n_nodes - 1
        if not (0 <= self.discharge_outlet < self.n_nodes):
            raise ValueError("discharge_outlet index out of range.")

# --------- Core model (1D multi-node) ---------

def simulate_stratified_tank(
    cfg: TankConfig,
    Tin_charge: np.ndarray,      # [K] charging inlet temperature time series (used when m_flow>0)
    m_flow: np.ndarray,          # [kg/s] signed mass flow: +charge, -discharge, 0 idle
    T_amb: np.ndarray,           # [K] ambient temperature time series for losses
    T_init: np.ndarray,          # [K] initial vertical profile (n_nodes,)
) -> Dict[str, np.ndarray]:
    """
    1D multi-node stratified sensible storage model (finite volumes).

    Flow convention:
      - m_flow[t] > 0 : CHARGING (mass enters at cfg.charge_inlet with Tin_charge[t])
      - m_flow[t] < 0 : DISCHARGING (mass leaves at cfg.discharge_outlet)
      - m_flow[t] = 0 : IDLE

    Returns:
      dict with keys:
        'T'          : (n_steps+1, n_nodes) temperature profile history
        'T_top'      : (n_steps+1,) top node temperature
        'T_bot'      : (n_steps+1,) bottom node temperature
        'Q_loss'     : (n_steps,) instantaneous loss power [W]
        'exec_time_s': float, execution time in seconds
    """
    t0 = time.perf_counter()

    # --- Validate config and inputs ---
    cfg.validate()
    n = cfg.n_steps

    Tin_charge = _ensure_1d(Tin_charge, "Tin_charge")
    m_flow     = _ensure_1d(m_flow, "m_flow")
    T_amb      = _ensure_1d(T_amb, "T_amb")
    T_init     = np.asarray(T_init, dtype=float)

    _check_len("Tin_charge", Tin_charge, n)
    _check_len("m_flow", m_flow, n)
    _check_len("T_amb", T_amb, n)

    if T_init.shape != (cfg.n_nodes,):
        raise ValueError(f"T_init must be shape ({cfg.n_nodes},), got {T_init.shape}.")

    # --- Geometry (constant properties) ---
    dz = cfg.height / cfg.n_nodes
    area_cross = cfg.volume / cfg.height      # A = V/H
    vol_node = area_cross * dz
    mass_node = cfg.rho * vol_node            # per node
    ua_node = cfg.ua_loss / cfg.n_nodes if cfg.ua_loss > 0 else 0.0
    G_cond = cfg.k_cond * area_cross / dz if cfg.k_cond > 0 else 0.0  # W/K per interface

    # --- State arrays ---
    T = np.empty((n + 1, cfg.n_nodes), dtype=float)
    T[0, :] = T_init
    Q_loss = np.zeros(n, dtype=float)

    # --- Time marching ---
    for t in range(n):
        Tn = T[t, :].copy()

        # 1) Side losses (explicit Euler)
        if ua_node > 0.0:
            dT_loss = (ua_node * (Tn - T_amb[t]) * cfg.dt) / (mass_node * cfg.cp)
            Tn = Tn - dT_loss
            Q_loss[t] = np.sum(ua_node * (T[t, :] - T_amb[t]))

        # 2) Axial conduction (explicit Euler)
        if G_cond > 0.0:
            dT_cond = np.zeros_like(Tn)
            # bottom
            dT_cond[0] = (G_cond * (Tn[1] - Tn[0]) * cfg.dt) / (mass_node * cfg.cp)
            # interior
            dT_cond[1:-1] = (G_cond * ((Tn[0:-2] - Tn[1:-1]) + (Tn[2:] - Tn[1:-1])) * cfg.dt) / (mass_node * cfg.cp)
            # top
            dT_cond[-1] = (G_cond * (Tn[-2] - Tn[-1]) * cfg.dt) / (mass_node * cfg.cp)
            Tn = Tn + dT_cond

        mf = _as_scalar(m_flow[t])

        # 3) Charging (+) / Discharging (-) / Idle (0)
        if mf > 0.0:
            # CHARGING at inlet node (instantaneous mixing)
            i_in = cfg.charge_inlet
            alpha = (mf * cfg.dt) / (mass_node + mf * cfg.dt)
            Tn[i_in] = (1 - alpha) * Tn[i_in] + alpha * _as_scalar(Tin_charge[t])

            # Upward advection proxy from inlet towards the top
            if i_in < cfg.n_nodes - 1:
                frac = min(1.0, (mf * cfg.dt) / mass_node)
                for i in range(i_in, cfg.n_nodes - 1):
                    Tn[i + 1] = (1 - frac) * Tn[i + 1] + frac * Tn[i]

        elif mf < 0.0:
            # DISCHARGING at outlet node; pull from flow direction side to preserve stratification
            i_out = cfg.discharge_outlet
            m_out = abs(mf)
            frac = min(1.0, (m_out * cfg.dt) / mass_node)

            if i_out > 0:
                # top outlet: pull from below upward
                Tn[i_out] = (1 - frac) * Tn[i_out] + frac * Tn[i_out - 1]
                for i in range(i_out - 1, 0, -1):
                    Tn[i] = (1 - frac) * Tn[i] + frac * Tn[i - 1]
            else:
                # bottom outlet: pull from above downward
                if cfg.n_nodes > 1:
                    Tn[i_out] = (1 - frac) * Tn[i_out] + frac * Tn[i_out + 1]
                    for i in range(0, cfg.n_nodes - 2):
                        Tn[i] = (1 - frac) * Tn[i] + frac * Tn[i + 1]

        else:
            # IDLE: no mass exchange
            pass

        # Optional physical clamp (disabled by default)
        # Tn = np.clip(Tn, 0.0, 473.15)

        T[t + 1, :] = Tn

    exec_time_s = time.perf_counter() - t0

    return {
        "T": T,
        "T_top": T[:, -1],
        "T_bot": T[:, 0],
        "Q_loss": Q_loss,
        "exec_time_s": exec_time_s,
    }
