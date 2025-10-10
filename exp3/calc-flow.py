# -*- coding: utf-8 -*-

import CoolProp.CoolProp as CP
from math import pi, e, sqrt
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

tube_diameter=74.7/1000               # Tubulation diameter [m]
tubulation_diameter=47.3/1000         # Orifice plate diameter [m]
pressure_0=692*133.332                # Air pressure [Pa]
temperature_0=24+273.15             # Air temperature [K]
work_fluid='air'                      # Work fluid (On tubulation)
metric_fluid='water'                  # Metric's fluid (On manometer)
g=9.78                                # Acceleration of Gravity [m/s^2]

betha = tubulation_diameter/tube_diameter                                           # Betha factor (No dimension)
tubulation_area = pi*pow(tubulation_diameter, 2)/4                                  # Area on orifice plate [m^2]
air_rho = CP.PropsSI("D", "T", temperature_0, "P", pressure_0, work_fluid)              # Air density [kg/m^3]
air_mu = CP.PropsSI("V", "T", temperature_0, "P", pressure_0, work_fluid)               # Air dynamic viscosity 
water_rho = CP.PropsSI("D","T",temperature_0,"P",pressure_0,'water')                    # Water density [kg/m^3]

data = {
    "30": [8, 20.1, 47, 53.5, 33.5],
    "60": [36, 76, 112, 158.5, 204.5]
} # Data: {Frequency: (Delta h for pressure difference, Delta h for static pressure) in [mmH20]}

for d in data: # Delta h in meter to delta P in Pa
    for num, p in enumerate(data[d]):
        data[d][num] *= water_rho*g/1000 # Conversion of water milimiter pressure to Pa

L1=1
L2 = 0.47
M2=2*L2/(1-betha)
results_x=[]
results_y=[]
results={}
for element in data:
    for num, p in enumerate(data[d]):
        dp=data[element][num]
        error=0.01        
        delta_re=0.02        
        Re=100          
        i=0
        while delta_re>error:
            A=pow(19e3*betha/Re, 0.8)
            C=0.5961 + 0.026*pow(betha, 2)-0.216*pow(betha, 8)+0.000521*pow(pow(10, 6)*betha/Re, 0.7) + (0.0188+0.0063*A)*pow(pow(10,6)/Re, 0.3)*pow(betha, 3.5)+(0.043+0.080*pow(e,-10*L1) - 0.0123*pow(e, -7*L1))*(1-0.11*A)*(pow(betha,4)/(1-pow(betha,4))) - 0.031*(M2 - 0.2*pow(M2, 1.1))*pow(betha, 1.3)
            massic_flow_rate=(C*tubulation_area*(2*air_rho*dp)**(1/2))/((1-betha**4)**(1/2))
            Re_new=4*massic_flow_rate/(pi*tube_diameter*air_mu)
            delta_re=abs(Re_new-Re)
            Re=Re_new
            i+=1
        volumetric_rate=massic_flow_rate*3600/air_rho             # Volumetric flow rate [m^3/h]
        if float(element) not in results:
            results[float(element)] = [volumetric_rate]
        else: results[float(element)].append(volumetric_rate)
print(results)