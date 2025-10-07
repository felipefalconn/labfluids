# -*- coding: utf-8 -*-
import numpy as np
from math import pi, e
import CoolProp.CoolProp as CP

# ---------- Configuração nominal (valores fornecidos) ----------
tube_diameter_nom = 74.7/1000          # m
tubulation_diameter_nom = 47.3/1000    # m
pressure_0_nom = 689 * 133.332         # Pa  (689 mmHg -> Pa)
temperature_0_nom = 25.2 + 273.15      # K
work_fluid = 'air'
g = 9.78                               # m/s^2

# Incertezas (desvios padrão)
sigma_tube_diameter = 0.05/1000        # m (0.05 mm)
sigma_tubulation_diameter = 0.05/1000  # m (0.05 mm)
sigma_pressure_0 = 0.5 * 133.332       # Pa (0.5 mmHg -> Pa)
sigma_temperature_0 = 0.5              # K  (0.5 °C)
sigma_dp_mm = 0.5                       # mm (incerteza em Delta h, antes da conversão)

# Dados originais (Δh em mm)
data = {
    "15": [12.5, 4.5],
    "30": [55.5, 14.5],
    "45": [123.5, 35.5],
    "60": [222.5, 65.5]
}

# Parâmetros do algoritmo iterativo (copiados do seu código)
L1 = 1
L2 = 0.47
results = {}

# Monte Carlo settings
n_samples = 5000   # aumentar para incerteza mais estável (ex.: 10000)
rng = np.random.default_rng(42)  # semente para reprodutibilidade

def compute_flow_rate(tube_diameter, tubulation_diameter, pressure_0, temperature_0, dp_mm):
    """Calcula vazão volumétrica (m^3/h) dado um conjunto de entradas (tudo em float)."""
    betha = tubulation_diameter / tube_diameter
    tubulation_area = pi * tubulation_diameter**2 / 4.0

    # densidades e viscosidade com CoolProp
    air_rho = CP.PropsSI("D", "T", temperature_0, "P", pressure_0, work_fluid)
    air_mu  = CP.PropsSI("V", "T", temperature_0, "P", pressure_0, work_fluid)
    water_rho = CP.PropsSI("D", "T", temperature_0, "P", pressure_0, "water")

    # converte dp (mm de coluna d'água) para Pa usando water_rho
    dp_pa = (dp_mm / 1000.0) * water_rho * g

    # iteração para Re (cópia do seu loop)
    M2 = 2 * L2 / (1 - betha)
    error = 0.01
    delta_re = 0.02
    Re = 100.0
    i = 0
    while delta_re > error and i < 300:
        A = (19e3 * betha / Re) ** 0.8
        C = (0.5961 + 0.026 * betha**2 - 0.216 * betha**8
             + 0.000521 * ((1e6 * betha / Re) ** 0.7)
             + (0.0188 + 0.0063 * A) * ((1e6 / Re) ** 0.3) * (betha ** 3.5)
             + (0.043 + 0.080 * e**(-10 * L1) - 0.0123 * e**(-7 * L1)) * (1 - 0.11 * A) * (betha**4 / (1 - betha**4))
             - 0.031 * (M2 - 0.2 * M2**1.1) * betha**1.3)
        massic_flow_rate = (C * tubulation_area * (2 * air_rho * dp_pa)**0.5) / ((1 - betha**4)**0.5)
        Re_new = 4 * massic_flow_rate / (pi * tube_diameter * air_mu)
        delta_re = abs(Re_new - Re)
        Re = Re_new
        i += 1

    volumetric_rate = massic_flow_rate * 3600.0 / air_rho  # m^3/h
    return volumetric_rate, massic_flow_rate, Re, i, C

# --- Calculo nominal (sem amostragem) para referência ---
for freq, vals in data.items():
    dp_mm_nom = vals[0]   # só usamos dp (Δh) que corresponde à pressão diferencial
    vol_nom, mass_nom, Re_nom, iters_nom, C_nom = compute_flow_rate(
        tube_diameter_nom,
        tubulation_diameter_nom,
        pressure_0_nom,
        temperature_0_nom,
        dp_mm_nom
    )
    results[freq] = {
        "nominal_m3h": vol_nom,
        "nominal_mass": mass_nom,
        "nominal_Re": Re_nom,
        "C_nominal": C_nom,
        "iters_nominal": iters_nom
    }

# --- Monte Carlo: amostragem das variáveis incertas e estimativa da incerteza das vazões ---
for freq, vals in data.items():
    dp_mm_nom = vals[0]   # mm
    samples_vol = np.empty(n_samples, dtype=float)

    # gerar amostras das variáveis incertas
    tube_samples = rng.normal(tube_diameter_nom, sigma_tube_diameter, size=n_samples)
    tubulation_samples = rng.normal(tubulation_diameter_nom, sigma_tubulation_diameter, size=n_samples)
    pressure_samples = rng.normal(pressure_0_nom, sigma_pressure_0, size=n_samples)
    temp_samples = rng.normal(temperature_0_nom, sigma_temperature_0, size=n_samples)
    dp_mm_samples = rng.normal(dp_mm_nom, sigma_dp_mm, size=n_samples)

    # loop de Monte Carlo
    for k in range(n_samples):
        try:
            vol_k, _, _, _, _ = compute_flow_rate(
                tube_samples[k],
                tubulation_samples[k],
                pressure_samples[k],
                temp_samples[k],
                dp_mm_samples[k]
            )
        except Exception as e:
            # caso algo falhe (valores inválidos), marca NaN e seguirá
            vol_k = np.nan
        samples_vol[k] = vol_k

    # descartar NaNs (se ocorreram)
    samples_vol = samples_vol[~np.isnan(samples_vol)]
    mean_vol = np.mean(samples_vol)
    std_vol = np.std(samples_vol, ddof=1)   # desvio padrão amostral
    rel_percent = (std_vol / mean_vol) * 100.0 if mean_vol != 0 else np.nan

    results[freq].update({
        "mc_mean_m3h": mean_vol,
        "mc_std_m3h": std_vol,
        "mc_rel_percent": rel_percent,
        "mc_samples_used": len(samples_vol)
    })

# --- Exibir resultados ---
print("Resultados por frequência (vazão volumétrica):")
print("Freq [Hz] | Nominal [m^3/h] | MC mean [m^3/h] ± std | incerteza rel (%) | samples_used")
for freq in sorted(results, key=lambda x: int(x)):
    r = results[freq]
    print(f"{freq:>6} | {r['nominal_m3h']:.6f}       | {r['mc_mean_m3h']:.6f} ± {r['mc_std_m3h']:.6f} | {r['mc_rel_percent']:.2f}% | {r['mc_samples_used']}")

# Opcional: acessar o dicionário `results` para salvar/plotar/exportar.
