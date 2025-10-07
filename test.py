import math
from CoolProp.CoolProp import PropsSI

# Parâmetros fixos
f = 6.49  # Calculado via fórmula de Barr (1931)
g = 9.78
m_e = 0.01945  # kg
D_e = 0.03025  # m (diâmetro da esfera)
D_vazao = 0.0338  # m (diâmetro para cálculo da vazão)

# Condições ambientais
pressao = 694 * 133.322  # Conversão de mmHg para Pa
temperatura = 25.1 + 273.15  # Conversão para Kelvin

# Cálculo das propriedades do ar usando CoolProp
rho_ar = PropsSI('D', 'P', pressao, 'T', temperatura, 'Air')  # Densidade [kg/m³]
mu_ar = PropsSI('V', 'P', pressao, 'T', temperatura, 'Air')  # Viscosidade [Pa·s]

# Geometria
A_e = math.pi * (D_e ** 2) / 4  # Área frontal da esfera [m²]
V_e = math.pi * (D_e ** 3) / 6  # Volume da esfera [m³]
A_vazao = math.pi * (D_vazao ** 2) / 4  # Área para vazão [m²]

# Função para cálculo do coeficiente de arrasto (Equação 8)
def calcular_Cd(Re):
    termo1 = 24 / Re
    termo2 = (2.6 * (Re / 5.0)) / (1 + (Re / 5.0) ** 1.52)
    termo3 = (0.411 * (Re / 2.63e5) ** -7.94) / (1 + (Re / 2.63e5) ** -8.00)
    termo4 = (0.25 * (Re / 1e6)) / (1 + Re / 1e6)
    return termo1 + termo2 + termo3 + termo4

# Processo iterativo
Re_old = 100
tolerancia = 1

for _ in range(20):  # Máximo de 20 iterações
    Cd = calcular_Cd(Re_old)
    Cdm = Cd * f**2
    
    # Correção: Toda a expressão dentro da raiz quadrada
    v = math.sqrt((2 * g * (m_e - rho_ar * V_e)) / (Cdm * rho_ar * A_e))
    Re_new = (rho_ar * v * D_e) / mu_ar
    
    if abs(Re_new - Re_old) < tolerancia:
        break
    
    Re_old = Re_new

vazao_massica = rho_ar * v * A_vazao

# Exibição dos resultados
print(f'Propriedades do ar:')
print(f'Densidade (): {rho_ar:.4f} kg/m³')
print(f'Viscosidade (): {mu_ar:.2e} Pa·s\n')
print(f'Resultados finais:')
print(f'Re final: {Re_new:.0f}')
print(f'Velocidade: {v:.2f} m/s')
print(f'Vazão mássica: {vazao_massica:.5f} kg/s')
print(f'A_vazao: {A_vazao}')