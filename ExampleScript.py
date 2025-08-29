# -*- coding: utf-8 -*-
"""
Demonstrates the 1D multi-node stratified tank with realistic operation modes:
charging, idle, discharging, and a combined schedule.

Developed on Fri Aug 29 18:25:43 2025

@author: Dr. Hakan İbrahim Tol (thanks to ChatGPT)
"""

"""
Four scenarios for the 1D multi-node stratified tank:
  (1) Full charging from a cool initial state (bottom inlet, hot water)
  (2) Discharging from a charged state (top outlet)
  (3) Idle from a stratified initial profile (losses + conduction)
  (4) Charging -> Idle -> Discharging with changing flow rates

Requires: storage1d.py in same folder.

Run directly.

"""

import numpy as np
import matplotlib.pyplot as plt
from storage1d import TankConfig, simulate_stratified_tank


# ----------------------------- Plot helpers -----------------------------

def build_time_axes(n_steps, dt):
    """Return (time_edges_min, time_steps_min)."""
    time_edges = np.arange(n_steps + 1) * (dt / 60.0)  # length n+1
    time_steps = np.arange(n_steps) * (dt / 60.0)      # length n
    return time_edges, time_steps

def plot_top_bottom(T_top, T_bot, time_edges, title):
    plt.figure()
    plt.plot(time_edges, T_top - 273.15, label="Top temperature")
    plt.plot(time_edges, T_bot - 273.15, label="Bottom temperature")
    plt.xlabel("Time [min]")
    plt.ylabel("Temperature [°C]")
    plt.title(title)
    plt.legend()
    plt.grid(True)

def plot_temp_field(T_full, cfg, time_edges, title):
    plt.figure()
    extent = [time_edges[0], time_edges[-1], 0.0, cfg.height]
    T_c = (T_full - 273.15)[:, ::-1]  # flip so bottom at y=0 in the image
    plt.imshow(
        T_c.T, aspect='auto', origin='lower', extent=extent, interpolation='bilinear'
    )
    plt.colorbar(label="Temperature [°C]")
    plt.xlabel("Time [min]")
    plt.ylabel("Tank height [m]")
    plt.title(title)

def plot_mflow(m_flow, time_steps, title):
    plt.figure()
    plt.step(time_steps, m_flow, where="post")
    plt.xlabel("Time [min]")
    plt.ylabel("Mass flow [kg/s] (+: charge, −: discharge)")
    plt.title(title)
    plt.grid(True)

def plot_profiles_height_on_x(T_full, cfg, times_min_selected, title):
    """
    x = height [m], y = temperature [°C]; snapshots at selected minutes.
    """
    dz = cfg.height / cfg.n_nodes
    z_centers = np.linspace(dz/2, cfg.height - dz/2, cfg.n_nodes)
    dt_min = cfg.dt / 60.0

    plt.figure()
    for tmin in times_min_selected:
        idx = int(round(tmin / dt_min))
        idx = max(0, min(idx, T_full.shape[0] - 1))
        T_profile_C = T_full[idx, :] - 273.15
        plt.plot(z_centers, T_profile_C, label=f"{tmin:.0f} min")

    plt.xlabel("Tank height [m]")
    plt.ylabel("Temperature [°C]")
    plt.title(title)
    plt.grid(True)
    plt.legend()


# ----------------------------- Scenarios -----------------------------

def scenario_full_charge(cfg):
    """
    (1) Full charging from cool initial state:
        - Uniformly cool tank
        - Positive constant m_flow (charge)
        - Hot inlet at bottom
    """
    n = cfg.n_steps
    Tin_charge = np.full(n, 60.0 + 273.15)   # 60°C
    T_amb      = np.full(n, 20.0 + 273.15)   # 20°C
    T_init     = np.full(cfg.n_nodes, 20.0 + 273.15)  # cool, uniform
    m_flow     = np.full(n, 0.10)            # kg/s charging

    out = simulate_stratified_tank(cfg, Tin_charge, m_flow, T_amb, T_init)
    return out, m_flow, "Scenario 1 — Full Charging (cool → hot, bottom inlet)"


def scenario_discharge_from_charged(cfg, T_init_from=None):
    """
    (2) Discharging from charged state:
        - If T_init_from provided (e.g., final state of full charge), use it.
        - Else create a hot-stratified initial profile.
        - Negative constant m_flow (discharge) at top outlet
    """
    n = cfg.n_steps
    Tin_charge = np.full(n, 60.0 + 273.15)   # not used when m_flow<0; still required
    T_amb      = np.full(n, 20.0 + 273.15)

    if T_init_from is None:
        # Hotter at top, cooler at bottom (e.g., 55→65°C)
        T_init = np.linspace(55.0, 65.0, cfg.n_nodes) + 273.15
    else:
        T_init = T_init_from.copy()

    m_flow = np.full(n, -0.08)               # kg/s discharging

    out = simulate_stratified_tank(cfg, Tin_charge, m_flow, T_amb, T_init)
    return out, m_flow, "Scenario 2 — Discharging from Charged State (top outlet)"


def scenario_idle_stratified(cfg):
    """
    (3) Idle from a stratified initial profile:
        - No mass flow (m_flow = 0)
        - Include losses and/or conduction effects if configured
    """
    n = cfg.n_steps
    Tin_charge = np.full(n, 60.0 + 273.15)  # unused in idle
    T_amb      = np.full(n, 20.0 + 273.15)

    # Stratified: e.g., 30→50°C from bottom to top
    T_init = np.linspace(30.0, 50.0, cfg.n_nodes) + 273.15
    m_flow = np.zeros(n)

    out = simulate_stratified_tank(cfg, Tin_charge, m_flow, T_amb, T_init)
    return out, m_flow, "Scenario 3 — Idle (no flow) from Stratified Profile"


def scenario_charge_idle_discharge(cfg):
    """
    (4) Charging → Idle → Discharging with changing flow rates:
        - 0–30 min: charge at 0.12 kg/s
        - 30–45 min: idle
        - 45–120 min: discharge at 0.06 kg/s
    """
    n = cfg.n_steps
    Tin_charge = np.full(n, 60.0 + 273.15)
    T_amb      = np.full(n, 20.0 + 273.15)
    T_init     = np.full(cfg.n_nodes, 22.0 + 273.15)  # slightly cool start

    m_flow = np.zeros(n)
    steps_30 = int(30 * 60 / cfg.dt)
    steps_45 = int(45 * 60 / cfg.dt)

    m_flow[:steps_30] = 0.12    # charge
    m_flow[steps_30:steps_45] = 0.0   # idle
    m_flow[steps_45:] = -0.06   # discharge

    out = simulate_stratified_tank(cfg, Tin_charge, m_flow, T_amb, T_init)
    return out, m_flow, "Scenario 4 — Charge → Idle → Discharge (varying flow)"


# ----------------------------- Main runner -----------------------------

def run_and_plot(cfg, out, m_flow, scenario_title):
    T = out["T"]
    T_top = out["T_top"]
    T_bot = out["T_bot"]
    exec_time_s = out["exec_time_s"]

    print(f"{scenario_title}: exec_time_s = {exec_time_s:.6f} s")
    print(f"  Final top T:    {T_top[-1]-273.15:.2f} °C")
    print(f"  Final bottom T: {T_bot[-1]-273.15:.2f} °C")

    time_edges, time_steps = build_time_axes(cfg.n_steps, cfg.dt)

    # Plots
    plot_top_bottom(T_top, T_bot, time_edges, scenario_title + " — Top/Bottom T")
    plot_temp_field(T, cfg, time_edges, scenario_title + " — Temperature Field")
    plot_mflow(m_flow, time_steps, scenario_title + " — Mass Flow Schedule")

    # Profiles: pick a few representative snapshots
    # (start, early, mid, late, end) — adjust per scenario duration
    total_min = time_edges[-1]
    snapshots = [0, 15, 30, 60, total_min]
    plot_profiles_height_on_x(T, cfg, snapshots, scenario_title + " — Profiles T(z)")


def main():
    # Global configuration (adjust as needed)
    cfg = TankConfig(
        n_nodes=20,
        height=2.0,
        volume=1.0,
        ua_loss=5.0,      # W/K total losses
        k_cond=0.10,      # W/m/K axial conduction
        dt=60.0,          # 60 s
        n_steps=120,      # 120 minutes total
        charge_inlet=0,   # bottom inlet
        discharge_outlet=-1   # top outlet
    )

    # (1) Full charging from cool initial state
    out1, m1, title1 = scenario_full_charge(cfg)
    run_and_plot(cfg, out1, m1, title1)

    # (2) Discharging from charged state (use final T of scenario 1)
    T_init2 = out1["T"][-1, :]
    out2, m2, title2 = scenario_discharge_from_charged(cfg, T_init_from=T_init2)
    run_and_plot(cfg, out2, m2, title2)

    # (3) Idle from stratified profile
    out3, m3, title3 = scenario_idle_stratified(cfg)
    run_and_plot(cfg, out3, m3, title3)

    # (4) Charge -> Idle -> Discharge with changing flow
    out4, m4, title4 = scenario_charge_idle_discharge(cfg)
    run_and_plot(cfg, out4, m4, title4)

    plt.show()


if __name__ == "__main__":
    main()
