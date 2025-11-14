from CoolProp.CoolProp import PropsSI
from uncertainties import ufloat
import numpy as np
import pandas as pd

P0_u = 0.5*0.133322
t_u = 0.5/100
P0 = ufloat(698*0.13322, P0_u)
T0 = ufloat(21+273.15, 0.5)
g = 9.78
rho_w = PropsSI("D","T",T0.n,"P",P0.n*1000,"water")

dPR_bar = 0.05
dPS_bar = 0.05
dPS_cm = 1.0

data_b1 = {
    "PR":[2,2.1,2.3,2.4,2.5],
    "PS":[-0.1,-0.1,-0.1,0,0],
    "Q":[9000,8500,7500,6500,5500],
    "m":[ufloat(0.395,t_u)]*5
}

data_b2 = {
    "PR":[1.9,1.95,2.15,2.3,2.4],
    "PS":[325,339,352,381,394],
    "Q":[9000,8500,7500,6500,5500],
    "m":[ufloat(0.48,t_u)]*5
}

data_series = {
    "PR":[2.9,3.3,3.9,4.4,4.7],
    "PS":[265,288,332,365,401],
    "Q":[10500,9500,8500,7500,6500],
    "m":[ufloat(1.493,t_u)]*5
}

data_parallel = {
    "PR1":[2.6,2.7,2.8,2.9,2.9],
    "PS1":[0,0,0,0,0],
    "PR2":[2.5,2.6,2.65,2.7,2.7],
    "PS2":[404,412,424,428,432],
    "Q":[10000,8000,7000,6000,5000],
    "m":[ufloat(0.085,t_u)]*5
}

def tabela_bomba1(data):
    Q = np.array(data["Q"])/3.6e6
    H_u = []
    for pr, ps, m in zip(data["PR"], data["PS"], data["m"]):
        PR_u = ufloat(pr, dPR_bar)*101300
        PS_u = ufloat(ps, dPS_bar)*101300
        H = PR_u/(rho_w*g) - PS_u/(rho_w*g) + m
        H_u.append(H)
    H_nom = [h.n for h in H_u]
    H_std = [h.s for h in H_u]
    return pd.DataFrame({"Q [m³/s]":Q,"H [m]":H_nom,"ΔH [m]":H_std})

def tabela_bomba2(data):
    Q = np.array(data["Q"])/3.6e6
    H_u = []
    for pr, ps_cm, m in zip(data["PR"], data["PS"], data["m"]):
        PR_u = ufloat(pr, dPR_bar)*101300
        PS_u = ufloat(ps_cm, dPS_cm)/100*rho_w*g
        H = PR_u/(rho_w*g) - PS_u/(rho_w*g) + m
        H_u.append(H)
    H_nom = [h.n for h in H_u]
    H_std = [h.s for h in H_u]
    return pd.DataFrame({"Q [m³/s]":Q,"H [m]":H_nom,"ΔH [m]":H_std})

def tabela_series(data):
    Q = np.array(data["Q"])/3.6e6
    H_u = []
    for pr, ps_cm, m in zip(data["PR"], data["PS"], data["m"]):
        PR_u = ufloat(pr, dPR_bar)*101300
        PS_u = ufloat(ps_cm, dPS_cm)/100*rho_w*g
        H = PR_u/(rho_w*g) - PS_u/(rho_w*g) + m
        H_u.append(H)
    H_nom = [h.n for h in H_u]
    H_std = [h.s for h in H_u]
    return pd.DataFrame({"Q [m³/s]":Q,"H [m]":H_nom,"ΔH [m]":H_std})

def tabela_parallel(data):
    Q = np.array(data["Q"])/3.6e6
    H_u = []
    for pr1, pr2, ps1_cm, ps2_cm, m in zip(
        data["PR1"], data["PR2"], data["PS1"], data["PS2"], data["m"]
    ):
        PR1_u = ufloat(pr1, dPR_bar)*101300
        PR2_u = ufloat(pr2, dPR_bar)*101300
        PR_med = (PR1_u + PR2_u)/2
        PS1_u = ufloat(ps1_cm, dPS_cm)/100*rho_w*g
        PS2_u = ufloat(ps2_cm, dPS_cm)/100*rho_w*g
        PS_med = (PS1_u + PS2_u)/2
        H = PR_med/(rho_w*g) - PS_med/(rho_w*g) + m
        H_u.append(H)
    H_nom = [h.n for h in H_u]
    H_std = [h.s for h in H_u]
    return pd.DataFrame({"Q [m³/s]":Q,"H [m]":H_nom,"ΔH [m]":H_std})

tab_b1 = tabela_bomba1(data_b1)
tab_b2 = tabela_bomba2(data_b2)
tab_series = tabela_series(data_series)
tab_parallel = tabela_parallel(data_parallel)

print("\nBomba 1:")
print(tab_b1.to_string(index=False))
print("\nBomba 2:")
print(tab_b2.to_string(index=False))
print("\nSérie (experimental):")
print(tab_series.to_string(index=False))
print("\nParalelo (experimental):")
print(tab_parallel.to_string(index=False))
