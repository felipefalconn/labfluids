from CoolProp.CoolProp import PropsSI
from uncertainties import ufloat
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({
    "figure.figsize": (8, 5),
    "axes.grid": True,
    "grid.linestyle": ":",
    "grid.alpha": 0.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#444444",
    "axes.labelcolor": "#222222",
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "legend.frameon": False,
    "legend.fontsize": 9
})

c1 = "#1f77b4"
c2 = "#2ca02c"
c3 = "#7f7f7f"
c4 = "#d62728"
c5 = "#9467bd"
c6 = "#8c564b"

P0_u = 0.5*0.133322
t_u = 0.5/100
P0 = ufloat(698*0.13322, P0_u)
T0 = ufloat(21+273.15, 0.5)
g = 9.78
rho_w = PropsSI("D","T",T0.n,"P",P0.n*1000,'water')

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

dtseries = pd.DataFrame(data_series)
dtseries["PR"] = dtseries["PR"]*101300
dtseries["PS"] = (dtseries["PS"]/100)*rho_w*g
dtseries["Q"] = dtseries["Q"]/3.6e6
dtseries["H"] = (dtseries["PR"]/(rho_w*g) - dtseries["PS"]/(rho_w*g)) + [m.n for m in dtseries["m"]]

dtparallel = pd.DataFrame(data_parallel)
dtparallel["PR"] = (dtparallel["PR1"]*101300 + dtparallel["PR2"]*101300)/2
dtparallel["PS"] = ((dtparallel["PS1"]+dtparallel["PS2"])/2/100)*rho_w*g
dtparallel["Q"] = dtparallel["Q"]/3.6e6
dtparallel["H"] = (dtparallel["PR"]/(rho_w*g) - dtparallel["PS"]/(rho_w*g)) + [m.n for m in dtparallel["m"]]

dt1 = pd.DataFrame(data_b1)
dt1["PR"] = dt1["PR"]*101300
dt1["PS"] = dt1["PS"]*101300
dt1["Q"] = dt1["Q"]/3.6e6
dt1["H"] = (dt1["PR"]/(rho_w*g) - dt1["PS"]/(rho_w*g)) + [m.n for m in dt1["m"]]

dt2 = pd.DataFrame(data_b2)
dt2["PR"] = dt2["PR"]*101300
dt2["PS"] = (dt2["PS"]/100)*rho_w*g
dt2["Q"] = dt2["Q"]/3.6e6
dt2["H"] = (dt2["PR"]/(rho_w*g) - dt2["PS"]/(rho_w*g)) + [m.n for m in dt2["m"]]

dt1_h = dt1.sort_values("H")
dt2_h = dt2.sort_values("H")
dt1_q = dt1.sort_values("Q")
dt2_q = dt2.sort_values("Q")

H_min = max(dt1_h["H"].min(), dt2_h["H"].min())
H_max = min(dt1_h["H"].max(), dt2_h["H"].max())

H_fine = np.linspace(H_min, H_max, 400)
Q1_fine = np.interp(H_fine, dt1_h["H"], dt1_h["Q"])
Q2_fine = np.interp(H_fine, dt2_h["H"], dt2_h["Q"])
Q_tot = Q1_fine + Q2_fine

ordem = np.argsort(Q_tot)
Q_tot = Q_tot[ordem]
H_fine = H_fine[ordem]

Q_cut = 0.0022
Q_max = Q_tot.max()
Q_par_sample = np.linspace(Q_cut, Q_max, 5)
mask_dom = Q_tot >= Q_cut
Q_dom = Q_tot[mask_dom]
H_dom = H_fine[mask_dom]
H_par_sample = np.interp(Q_par_sample, Q_dom, H_dom)

fig, ax = plt.subplots()

ax.plot(dt1_q["Q"], dt1_q["H"], marker="o", linestyle="-", color=c1, label="Bomba 1")
ax.plot(dt2_q["Q"], dt2_q["H"], marker="s", linestyle="--", color=c2, label="Bomba 2")
ax.plot(dt1_q["Q"], dt1_q["H"]+dt2_q["H"], marker="D", linestyle="-", color=c3, label="Série (calc)")
ax.plot(Q_par_sample, H_par_sample, marker="o", linestyle="--", color=c4, label="Paralelo (calc, 5 pts)")
ax.plot(dtseries["Q"], dtseries["H"], marker="^", linestyle="-.", color=c5, label="Série (exp)")
ax.plot(dtparallel["Q"], dtparallel["H"], marker="+", linestyle="-", color=c6, label="Paralelo (exp)")

ax.set_title("Curvas H×Q das Bombas")
ax.set_xlabel("Vazão Q [m³/s]")
ax.set_ylabel("Carga H [m]")
ax.legend()
plt.tight_layout()
plt.show()