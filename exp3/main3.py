import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import CoolProp.CoolProp as CP
import numpy as np

P0 = 692 * 133.322
T0 = 24 + 273.15
work_fluid = 'air'
water_rho = CP.PropsSI('D', 'T', T0, 'P', P0, 'water')
g = 9.78
data = { # Frequency: [[Q, H]]
    '30': [
        [0.013810316398790047, 74.5],
        [0.021702305560092602, 71],
        [0.03298646710987445, 55.5],
        [0.03516572996249611, 52.5],
        [0.02791090564504904, 60]
    ],
    '60': [
        [0.028919412330791098, 297.5],
        [0.04182921309143935, 278],
        [0.05067775037053948, 259],
        [0.06018997743896055, 236],
        [0.06829421441945342, 214]
    ]
}
rateux = {'60': []}

for e in data:
    for num in range(0, len(data[e])):
        data[e][num][1] *= 1/1000  # mmH20 to mH20
        if int(e) == 30:
            N_calc = water_rho * g * data[e][num][0] * data[e][num][1] / 0.283
            data[e][num].append(N_calc)
            
            Qnew = 2 * data[e][num][0]
            Hnew = pow(2, 2) * data[e][num][1]
            Nnew = pow(2, 3) * N_calc
            rateux['60'].append([Qnew, Hnew, Nnew])
        else:
            data[e][num].append(water_rho * g * data[e][num][0] * data[e][num][1] / 0.311)

df_30 = pd.DataFrame(data['30'], columns=['Q [m^3/s]', 'H [m]', 'N [W]'])
df_60 = pd.DataFrame(data['60'], columns=['Q [m^3/s]', 'H [m]', 'N [W]'])
df_60rateux = pd.DataFrame(rateux['60'], columns=['Q [m^3/s]', 'H [m]', 'N [W]'])

df_30.sort_values(by='Q [m^3/s]', ascending=True, inplace=True)
df_60.sort_values(by='Q [m^3/s]', ascending=True, inplace=True)
df_60rateux.sort_values(by='Q [m^3/s]', ascending=True, inplace=True)

def suavizar_curvas(df):
    """Aplica interpolação Cubic Spline para suavizar H e N em relação a Q."""
    Q_orig = df['Q [m^3/s]'].values
    H_orig = df['H [m]'].values
    N_orig = df['N [W]'].values

    Q_suavizado = np.linspace(Q_orig.min(), Q_orig.max(), 100)

    H_spline = CubicSpline(Q_orig, H_orig)
    N_spline = CubicSpline(Q_orig, N_orig)

    H_suavizado = H_spline(Q_suavizado)
    N_suavizado = N_spline(Q_suavizado)
    
    return Q_suavizado, H_suavizado, N_suavizado

Q30s, H30s, N30s = suavizar_curvas(df_30)
Q60s, H60s, N60s = suavizar_curvas(df_60)
Qrateuxs, Hrateuxs, Nrateuxs = suavizar_curvas(df_60rateux)


plt.style.use('seaborn-v0_8-whitegrid') 
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

colors = ['#C0392B', '#2980B9', '#27AE60', '#f556e5']
markers = ['o', 's', 'D']
linestyles = ['-', '--', '-.']
marker_size = 5

ax1 = axes[0]
ax1.plot(Q30s, H30s, linestyle=linestyles[0], color=colors[0], label='30 Hz', linewidth=2.5)
ax1.plot(df_30['Q [m^3/s]'], df_30['H [m]'], 'o', color=colors[0], markersize=marker_size) 

ax1.plot(Q60s, H60s, linestyle=linestyles[0], color=colors[1], label='60 Hz', linewidth=2.5)
ax1.plot(df_60['Q [m^3/s]'], df_60['H [m]'], 's', color=colors[1], markersize=marker_size) 

ax1.plot(Qrateuxs, Hrateuxs, linestyle=linestyles[1], color=colors[3], label='60 Hz - Rateux', linewidth=2.5, alpha=0.7)
ax1.plot(df_60rateux['Q [m^3/s]'], df_60rateux['H [m]'], 'D', color=colors[3], markersize=marker_size, alpha=0.7) 

ax1.set_title('Curva de Desempenho: Carga vs. Vazão', fontsize=14, fontweight='bold')
ax1.set_xlabel('Vazão $Q$ [$m^3/s$]', fontsize=12)
ax1.set_ylabel('Altura Manométrica $H$ [$m$]', fontsize=12)
ax1.tick_params(axis='both', which='major', labelsize=10)
ax1.legend(loc='best', fontsize=10, frameon=True, shadow=True)


ax2 = axes[1]
ax2.plot(Q30s, N30s, linestyle=linestyles[0], color=colors[0], label='30 Hz', linewidth=2.5)
ax2.plot(df_30['Q [m^3/s]'], df_30['N [W]'], 'o', color=colors[0], markersize=marker_size)

ax2.plot(Q60s, N60s, linestyle=linestyles[0], color=colors[1], label='60 Hz', linewidth=2.5)
ax2.plot(df_60['Q [m^3/s]'], df_60['N [W]'], 's', color=colors[1], markersize=marker_size)

ax2.plot(Qrateuxs, Nrateuxs, linestyle=linestyles[1], color=colors[3], label='60 Hz - Rateux', linewidth=2.5, alpha=0.7)
ax2.plot(df_60rateux['Q [m^3/s]'], df_60rateux['N [W]'], 'D', color=colors[3], markersize=marker_size, alpha=0.7)

ax2.set_title('Curva de Desempenho: Potência vs. Vazão', fontsize=14, fontweight='bold')
ax2.set_xlabel('Vazão $Q$ [$m^3/s$]', fontsize=12)
ax2.set_ylabel('Potência $N$ [$W$]', fontsize=12)
ax2.tick_params(axis='both', which='major', labelsize=10)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax2.legend(loc='best', fontsize=10, frameon=True, shadow=True)

plt.tight_layout(pad=3.0)
plt.show()