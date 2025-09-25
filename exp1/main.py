# -*- coding: utf-8 -*-

import CoolProp.CoolProp as CP
from math import pi, e

tube_diameter=74.7/1000               # Tubulation diameter [m]
tubulation_diameter=47.3/1000         # Orifice plate diameter [m]
pressure_0=689*133.332                # Air pressure [Pa]
temperature_0=25.2+273.15             # Air temperature [K]
work_fluid='air'                      # Work fluid (On tubulation)
metric_fluid='water'                  # Metric's fluid (On manometer)
g=9.78                                # Acceleration of Gravity [m/s^2]

betha = tubulation_diameter/tube_diameter                                           # Betha factor (No dimension)
tubulation_area = pi*pow(tubulation_diameter, 2)/4                                  # Area on orifice plate [m^2]
air_rho = CP.PropsSI("D", "T", temperature_0, "P", pressure_0, work_fluid)              # Air density [kg/m^3]
air_mu = CP.PropsSI("V", "T", temperature_0, "P", pressure_0, work_fluid)               # Air dynamic viscosity 
water_rho = CP.PropsSI("D","T",temperature_0,"P",pressure_0,'water')                    # Water density [kg/m^3]

data = {
    "15": [12.5, 4.5],
    "30": [55.5, 14.5],
    "45": [123.5, 35.5],
    "60": [222.5, 65.5]
} # Data: {Frequency: (Delta h for pressure difference, Delta h for static pressure) in [m]}

for d in data: # Delta h in meter to delta P in Pa
    data[d][0]*=water_rho*g/1000 # Conversion of water milimiter pressure to Pa
    data[d][1]*=water_rho*g/1000 # Conversion of water milimiter pressure to Pa

L1=1
L2 = 0.47
M2=2*L2/(1-betha)
for element in data:
    dp=data[element][0]
    error=0.01        
    delta_re=0.02        
    Re=100          
    i=0
    while delta_re>error:
        A=pow(19e3*betha/Re, 0.8)
        C=0.5961 + 0.026*pow(betha, 2)-0.216*pow(betha, 8)+0.000521*pow(pow(10, 6)*betha/Re, 0.7) + (0.0188+0.0063*A)*pow(pow(10,6)/Re, 0.3)*pow(betha, 3.5)+(0.043+0.080*pow(e,-10*L1) - 0.0123*pow(e, -7*L1))*(1-0.11*A)*(pow(betha,4)/(1-pow(betha,4))) - 0.031*(M2 - 0.2*pow(M2, 1.1))*pow(betha, 1.3)
        m=(C*tubulation_area *(2*air_rho*dp)**(1/2))/((1-betha**4)**(1/2))
        Re_new=4*m/(pi*tube_diameter*air_mu)
        delta_re=abs(Re_new-Re)
        Re=Re_new
        i+=1
    volumetric_rate=m*3600/air_rho             # Volumetric flow rate [m^3/h]
    data[element].append({"Results": {"Volumetric Flow Rate": volumetric_rate, "Massic Flow Rate": m, "Re": Re, "Iterations": i, "C": C, "Betha": betha}})

for e in data:
    print(f'Frequency: {e}Hz — Data: {data[e]}')