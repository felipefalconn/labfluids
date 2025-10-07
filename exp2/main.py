# -*- coding: utf-8 -*-

from CoolProp.CoolProp import  PropsSI
from math import pi

d = 37/1000
D=10/1000           #diametro da esfera
m=5/1000           #massa da esfera
P=100E2             #pressao
T=25+273.15         #temperatura do ar

FR='air'            #fluido de trabalho
g=9.81              #aceleracao da gravidade

Vol=(4/3)*pi*(D/2)**3                       #volume da esfera
rho_esf=m/Vol                               #densidade da esfera
A=(pi*D**2)/4                             #area frontal
rho=PropsSI("D","T",T,"P",P,FR)             #densidade do ar
mu=PropsSI("V","T",T,"P",P,FR)             #densidade do ar
Fe=rho*g*Vol                               #forca de empuxo
Fp=rho_esf*g*Vol                            #forca peso da esfera
Fa=Fp-Fe                                    #forca de arrasto

erro=0.01        #criterio de convergencia para termino da iteracao
conv=100        #inicializacao do parametro utilizado para verificar a convergencia da iteracao.
V=100          #inicializacao da estimativa inicial
iter=0          #contador de iteracoes
f_correction = 1/((1-(d/D))*pow((1-0.5*pow(d/D, 2)), 0.5))
while conv>erro:
    Re=rho*V*D/mu
    C=(24/Re)+((2.6*(Re/5))/(1+(Re/5)**1.52))+((0.411*(Re/263000)**(-7.94))/(1+(Re/263000)**(-8)))+((0.25*(Re/1000000))/(1+(Re/1000000)))
    V_new=(2*Fa/(rho*C*A))**(1/2)     
    conv=abs(V_new-V)
    V=V_new
    iter=iter+1
Cdm = C*pow(f_correction, 2)
print("Velocidade=","%.3f" % V, "m/s")


       