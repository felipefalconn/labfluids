# -*- coding: utf-8 -*-

from CoolProp.CoolProp import  PropsSI
from math import pi, sqrt

tube_diameter=34.12/1000
ball_diameter=31.2/1000           
ball_mass=19.39/1000           
pressure=767.5*133.322          
temperature=22.2+273.15         

FR='air'            
g=9.78 

Vol=(4/3)*pi*pow(ball_diameter/2,3)                      
rho_ball=ball_mass/Vol                               
ball_area=pi*pow(ball_diameter, 2)/4                             
air_rho=PropsSI("D","T",temperature,"P",pressure,FR)             
mu=PropsSI("V","T",temperature,"P",pressure,FR)             
Fe=air_rho*g*Vol                         
Fp=rho_ball*g*Vol                            
Fa=Fp-Fe                                    

erro=0.01        
conv=100        
V=100          
iter=0
gamma = ball_diameter/tube_diameter
f_correction = 1/((1-pow(gamma, 2))*pow((1-0.5*pow(gamma, 2)), 0.5))
while conv>erro:
    Re=air_rho*V*ball_diameter/mu
    C=(24/Re)+((2.6*(Re/5))/(1+(Re/5)**1.52))+((0.411*(Re/263000)**(-7.94))/(1+(Re/263000)**(-8)))+((0.25*(Re/1000000))/(1+(Re/1000000)))
    Cdm = C*pow(f_correction, 2)
    V_new=sqrt(2*Fa/(air_rho*Cdm*ball_area))   
    conv=abs(V_new-V)
    V=V_new
    iter+=1
massic_flow = air_rho*V*pi*pow(tube_diameter, 2)*3600/4
print(f"Rho={air_rho}")
print(f"Air viscosity={mu}")
print(f"Re={Re} m/s")
print(f"Speed={V:.3f} m/s")
print(f"Massic flow rate {massic_flow} kg/h")