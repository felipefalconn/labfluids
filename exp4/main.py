from CoolProp.CoolProp import  PropsSI
from math import pi
from uncertainties import ufloat
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

P0_u = 0.5*0.133322
R_U = 0.5 /1000             
t_u = 0.5/100                
P0 = ufloat(698*0.13322, P0_u)
T0_u = 0.5
T0 = ufloat(21+273.15, T0_u)
g = 9.81
rho_w=PropsSI("D","T",T0.n,"P",P0.n*1000,'water')   

data_b1 = {
    "PR": [2,2.1,2.3, 2.4, 2.5], 
    "PS": [-0.1, -0.1, -0.1, 0, 0], 
    "Q": [9000, 8500, 7500, 6500, 5500], 
    "m": [ufloat(39.5/100, t_u), ufloat(39.5/100, t_u), ufloat(39.5/100, t_u), ufloat(39.5/100, t_u), ufloat(39.5/100, t_u)]}

data_b2 = {
    "PR": [1.9,1.95,2.15, 2.3, 2.4], 
    "PS": [325, 339, 352, 381, 394], 
    "Q": [9000, 8500, 7500, 6500, 5500], 
    "m": [ufloat(48/100, t_u), ufloat(48/100, t_u), ufloat(48/100, t_u), ufloat(48/100, t_u), ufloat(48/100, t_u)]}

manufacturer = {
    "Q_m3h": [8.6, 8.0, 7.4, 6.6, 6.0, 4.1],   
    "H_m":   [22.0, 23.0, 24.0, 25.0, 26.0, 28.0]
}

manufacturer["Q_m3s"] = [q / 3.6e3 for q in manufacturer["Q_m3h"]]

dt1=pd.DataFrame(data_b1)
dt1['PR'] = dt1['PR']*101300
dt1['PS'] = dt1['PS']*101300
dt1['Q'] = dt1['Q']/3.6e6
dt1['H'] = ((dt1['PR']/(rho_w*g) - dt1['PS']/(rho_w*g))) + [m.n for m in dt1['m']]

dt2=pd.DataFrame(data_b2)
dt2['PR'] = dt2['PR']*101300
dt2['PS'] = (dt2['PS']/100)*rho_w*g
dt2['Q'] = dt2['Q'] / 3.6e6
dt2['H'] = ((dt2['PR']/(rho_w*g) - dt2['PS']/(rho_w*g))) + [m.n for m in dt1['m']]

dt1_sorted = dt1.sort_values(by='Q')
dt2_sorted = dt2.sort_values(by='Q')

print(dt1_sorted)
print(dt2_sorted)

plt.figure(figsize=(10, 6))
plt.plot(dt1_sorted['Q'], dt1_sorted['H'], marker='s', linestyle='-', label='Bomba 1', color='black')
plt.plot(dt2_sorted['Q'], dt2_sorted['H'], marker='o', linestyle='--', label='Bomba 2', color='purple')
plt.plot(manufacturer["Q_m3s"], manufacturer["H_m"], marker='*', linestyle='-', color='blue', label='Bomba - Fabricante')
plt.title('Curvas HxQ das Bombas')
plt.xlabel('Vazão (Q) [m³/s]')
plt.ylabel('Carga Manométrica (H) [m]')
plt.legend()
plt.grid(True)
plt.show()