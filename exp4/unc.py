from CoolProp.CoolProp import PropsSI
from math import pi
from uncertainties import ufloat
import uncertainties.unumpy as unp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# -----------------------------
# Dados ambientais
# -----------------------------
P0_u = 0.5 * 0.133322        # incerteza da pressão ambiente (em atm)
R_U = 0.5 / 1000             # erro da régua [m]
t_u = 0.5 / 100              # erro da trena [m]

P0 = ufloat(698 * 0.13322, P0_u)
T0 = ufloat(21 + 273.15, 0.5)

g = 9.81
rho_w = PropsSI("D", "T", T0.n, "P", P0.n * 1000, "water")  # densidade da água [kg/m³]

# -----------------------------
# Dados experimentais
# -----------------------------
# bomba 1: PR, PS em atm; Q em L/h; m em m
data_b1 = {
    "PR": [2, 2.1, 2.3, 2.4, 2.5],
    "PS": [-0.1, -0.1, -0.1, 0, 0],
    "Q":  [9000, 8500, 7500, 6500, 5500],
    "m":  [ufloat(39.5/100, t_u)] * 5
}

# bomba 2: PR em atm; PS em cm H2O; Q em L/h; m em m
data_b2 = {
    "PR": [1.9, 1.95, 2.15, 2.3, 2.4],
    "PS": [325, 339, 352, 381, 394],  # cm H2O
    "Q":  [9000, 8500, 7500, 6500, 5500],
    "m":  [ufloat(48/100, t_u)] * 5
}

# incertezas das leituras
u_atm   = 0.05       # ±0,05 atm nos manômetros
u_cmH2O = 0.5        # ±0,5 cm na coluna de água

# -----------------------------
# Bomba 1 – com incertezas
# -----------------------------
PR1_atm = unp.uarray(data_b1["PR"], [u_atm]*5)    # atm
PS1_atm = unp.uarray(data_b1["PS"], [u_atm]*5)    # atm

PR1 = PR1_atm * 101300        # Pa
PS1 = PS1_atm * 101300        # Pa

Q1 = np.array(data_b1["Q"]) / 3.6e6  # m³/s

# transformar lista de ufloats numa uarray
m1 = unp.uarray(
    [m.n for m in data_b1["m"]],
    [m.s for m in data_b1["m"]]
)

H1 = (PR1/(rho_w*g) - PS1/(rho_w*g)) + m1  # H com incerteza

dt1 = pd.DataFrame({
    "Q":     Q1,
    "H_nom": unp.nominal_values(H1),
    "H_std": unp.std_devs(H1)
})

# -----------------------------
# Bomba 2 – com incertezas
# -----------------------------
PR2_atm = unp.uarray(data_b2["PR"], [u_atm]*5)
PR2 = PR2_atm * 101300             # Pa

PS2_cm = unp.uarray(data_b2["PS"], [u_cmH2O]*5)  # cm H2O
PS2 = (PS2_cm / 100.0) * rho_w * g               # m H2O → Pa

Q2 = np.array(data_b2["Q"]) / 3.6e6

m2 = unp.uarray(
    [m.n for m in data_b2["m"]],
    [m.s for m in data_b2["m"]]
)

H2 = (PR2/(rho_w*g) - PS2/(rho_w*g)) + m2

dt2 = pd.DataFrame({
    "Q":     Q2,
    "H_nom": unp.nominal_values(H2),
    "H_std": unp.std_devs(H2)
})

# ordenar pelas vazões
dt1_sorted = dt1.sort_values(by="Q")
dt2_sorted = dt2.sort_values(by="Q")

print("--- Bomba 1 ---")
print(dt1_sorted)
print("\n--- Bomba 2 ---")
print(dt2_sorted)

# -----------------------------
# Plot usando valor nominal
# -----------------------------
manufacturer = {
    "Q_m3h": [8.6, 8.0, 7.4, 6.6, 6.0, 4.1],
    "H_m":   [22.0, 23.0, 24.0, 25.0, 26.0, 28.0]
}
manufacturer["Q_m3s"] = [q / 3.6e3 for q in manufacturer["Q_m3h"]]

plt.figure(figsize=(10, 6))
plt.errorbar(dt1_sorted["Q"], dt1_sorted["H_nom"], yerr=dt1_sorted["H_std"],
             fmt="o-", label="Bomba 1")
plt.errorbar(dt2_sorted["Q"], dt2_sorted["H_nom"], yerr=dt2_sorted["H_std"],
             fmt="s--", label="Bomba 2")
plt.plot(manufacturer["Q_m3s"], manufacturer["H_m"],
         marker="^", linestyle="-", label="Bomba - Fabricante")

plt.title("Curvas H×Q das Bombas (com incerteza em H)")
plt.xlabel("Vazão Q [m³/s]")
plt.ylabel("Carga manométrica H [m]")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
