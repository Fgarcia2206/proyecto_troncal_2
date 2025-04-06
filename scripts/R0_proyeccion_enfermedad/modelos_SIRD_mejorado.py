import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from ipywidgets import interact, FloatSlider

# 1. Configuración inicial de datos
data = {
    'Ciudad': ['Norilsk', 'Kayerkan', 'Dudinka'],
    'Población': [180000, 25000, 20000],
    'Semanas': {
        1: {'S': [200000, 80000, 60000], 'I': [300, 10, 20], 'D': [0, 0, 0]},
        2: {'I': [1500, 700, 1200], 'D': [300, 120, 60]},
        3: {'I': [5000, 2300, 3000], 'D': [650, 1000, 900]}
    }
}

# 2. Modelo SIRD mejorado
def modelo_sird(y, t, beta, gamma, mu, N):
    S, I, R, D = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - (gamma + mu) * I
    dRdt = gamma * I
    dDdt = mu * I
    return [dSdt, dIdt, dRdt, dDdt]

# 3. Función para calibrar parámetros automáticamente
def calibrar_modelo(ciudad_idx):
    # Datos observados
    obs_I = [data['Semanas'][1]['I'][ciudad_idx], 
             data['Semanas'][2]['I'][ciudad_idx], 
             data['Semanas'][3]['I'][ciudad_idx]]
    obs_D = [0, 
             data['Semanas'][2]['D'][ciudad_idx], 
             data['Semanas'][3]['D'][ciudad_idx]]
    
    # Estimación inicial de CFR
    cfr = (obs_D[2] - obs_D[1]) / (obs_I[2] - obs_I[1])
    
    # Búsqueda de mejores parámetros (simplificado)
    best_R0 = 2.5  # Valor inicial (podría implementarse optimización)
    duracion_infeccion = 12
    mu = cfr / duracion_infeccion
    gamma = (1 - cfr) / duracion_infeccion
    beta = best_R0 * (gamma + mu)
    
    return beta, gamma, mu, cfr

# 4. Análisis de sensibilidad interactivo
@interact(
    ciudad=['Norilsk', 'Kayerkan', 'Dudinka'],
    R0=FloatSlider(min=1.0, max=5.0, step=0.1, value=2.5),
    duracion_infeccion=FloatSlider(min=5, max=20, step=1, value=12),
    CFR=FloatSlider(min=0.1, max=0.8, step=0.05, value=0.45)
)
def analisis_sensibilidad(ciudad, R0, duracion_infeccion, CFR):
    idx = data['Ciudad'].index(ciudad)
    N = data['Población'][idx]
    S0 = data['Semanas'][1]['S'][idx]
    I0 = data['Semanas'][1]['I'][idx]
    D0 = 0
    
    # Calcular parámetros
    mu = CFR / duracion_infeccion
    gamma = (1 - CFR) / duracion_infeccion
    beta = R0 * (gamma + mu)
    
    # Simular modelo
    t = np.linspace(0, 70, 70)
    solucion = odeint(modelo_sird, [S0, I0, 0, D0], t, args=(beta, gamma, mu, N))
    S, I, R, D = solucion.T
    
    # Puntos observados
    obs_t = [7, 14, 21]
    obs_I = [data['Semanas'][1]['I'][idx], 
             data['Semanas'][2]['I'][idx], 
             data['Semanas'][3]['I'][idx]]
    obs_D = [0, 
             data['Semanas'][2]['D'][idx], 
             data['Semanas'][3]['D'][idx]]
    
    # Gráfico
    plt.figure(figsize=(12, 6))
    plt.plot(t, S, 'b-', label='Susceptibles')
    plt.plot(t, I, 'r--', label='Infectados')
    plt.plot(t, R, 'g:', label='Recuperados')
    plt.plot(t, D, 'k-', label='Muertos')
    plt.scatter(obs_t, obs_I, color='r', marker='o', label='Infectados observados')
    plt.scatter(obs_t, obs_D, color='k', marker='x', label='Muertos observados')
    
    plt.title(f"Modelo SIRD - {ciudad}\nR0={R0:.2f}, Duración infección={duracion_infeccion}d, CFR={CFR:.0%}")
    plt.xlabel('Días desde inicio')
    plt.ylabel('Número de personas')
    plt.legend()
    plt.grid(True)
    plt.show()

# 5. Capacidad hospitalaria (ejemplo para Norilsk)
def modelo_sird_hospital(y, t, beta, gamma, mu, N, capacidad_hospital):
    S, I, R, D = y
    # Aumentar mortalidad si se excede la capacidad hospitalaria
    factor_sobrecarga = max(0, I - capacidad_hospital) / capacidad_hospital
    mu_efectivo = mu * (1 + factor_sobrecarga * 3)  # Mortalidad aumenta hasta 4x
    
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - (gamma + mu_efectivo) * I
    dRdt = gamma * I
    dDdt = mu_efectivo * I
    return [dSdt, dIdt, dRdt, dDdt]

# 6. Movilidad entre ciudades (modelo simplificado)
def modelo_sird_conectado(y, t, parametros):
    S1, I1, R1, D1, S2, I2, R2, D2 = y
    beta1, gamma1, mu1, N1, beta2, gamma2, mu2, N2, tasa_movilidad = parametros
    
    # Flujo entre ciudades
    flujo = tasa_movilidad * (I1/N1 - I2/N2)
    
    dS1dt = -beta1 * S1 * I1 / N1 - flujo * S1
    dI1dt = beta1 * S1 * I1 / N1 - (gamma1 + mu1) * I1 - flujo * I1
    dR1dt = gamma1 * I1 - flujo * R1
    dD1dt = mu1 * I1
    
    dS2dt = -beta2 * S2 * I2 / N2 + flujo * S1
    dI2dt = beta2 * S2 * I2 / N2 - (gamma2 + mu2) * I2 + flujo * I1
    dR2dt = gamma2 * I2 + flujo * R1
    dD2dt = mu2 * I2
    
    return [dS1dt, dI1dt, dR1dt, dD1dt, dS2dt, dI2dt, dR2dt, dD2dt]

# 7. Ejecutar análisis para cada ciudad
print("=== Análisis Individual por Ciudad ===")
for i, ciudad in enumerate(data['Ciudad']):
    beta, gamma, mu, cfr = calibrar_modelo(i)
    N = data['Población'][i]
    S0 = data['Semanas'][1]['S'][i]
    I0 = data['Semanas'][1]['I'][i]
    
    t = np.linspace(0, 70, 70)
    solucion = odeint(modelo_sird, [S0, I0, 0, 0], t, args=(beta, gamma, mu, N))
    S, I, R, D = solucion.T
    
    # Proyecciones
    casos_semana10 = int(I[-1] + R[-1] + D[-1])
    muertes_semana10 = int(D[-1])
    
    print(f"\n**{ciudad}**")
    print(f"- R0 estimado: {beta/(gamma + mu):.2f}")
    print(f"- CFR observado: {cfr:.1%}")
    print(f"- Proyección semana 10:")
    print(f"  Total infectados: {casos_semana10:,} ({casos_semana10/N:.1%} población)")
    print(f"  Muertes estimadas: {muertes_semana10:,} ({muertes_semana10/N:.1%} población)")