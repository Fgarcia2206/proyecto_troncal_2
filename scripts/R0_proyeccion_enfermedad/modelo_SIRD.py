import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Datos actualizados
data = {
    'Ciudad': ['Norilsk', 'Kayerkan', 'Dudinka'],
    'Población': [180000, 25000, 20000],
    'S0': [200000, 80000, 60000],  # Susceptibles iniciales (ajustados)
    'I0': [300, 10, 20],           # Infectados iniciales (Semana 1)
    'D0': [0, 0, 0],               # Muertos iniciales (Semana 1)
    'R0': [0, 0, 0]                # Recuperados iniciales (Semana 1)
}

df = pd.DataFrame(data)

# Ajustar valores inconsistentes (S0 + I0 no puede ser > Población)
for i in range(len(df)):
    if df.loc[i, 'S0'] + df.loc[i, 'I0'] > df.loc[i, 'Población']:
        df.loc[i, 'S0'] = df.loc[i, 'Población'] - df.loc[i, 'I0']

# Tasas observadas de mortalidad (CFR) de tus datos
cfr_norilsk = (650 - 300) / (5000 - 1500)  # (Muertos Sem3 - Sem2)/(Infectados Sem3 - Sem2)
cfr_kayerkan = (1000 - 120) / (2300 - 700)
cfr_dudinka = (900 - 60) / (3000 - 1200)
df['CFR'] = [cfr_norilsk, cfr_kayerkan, cfr_dudinka]

# Función del modelo SIRD
def modelo_sird(y, t, beta, gamma, mu, N):
    S, I, R, D = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - (gamma + mu) * I
    dRdt = gamma * I
    dDdt = mu * I
    return [dSdt, dIdt, dRdt, dDdt]

# Parámetros para cada ciudad
plt.figure(figsize=(15, 10))
colores = ['b', 'r', 'g']

for i in range(len(df)):
    # Parámetros específicos
    N = df.loc[i, 'Población']
    S0 = df.loc[i, 'S0']
    I0 = df.loc[i, 'I0']
    R0 = df.loc[i, 'R0']
    D0 = df.loc[i, 'D0']
    cfr = df.loc[i, 'CFR']
    
    # Estimación de parámetros
    duracion_infeccion = 12  # días (basado en la progresión de la enfermedad)
    gamma = (1 - cfr) / duracion_infeccion  # Tasa de recuperación
    mu = cfr / duracion_infeccion           # Tasa de mortalidad
    
    # Calcular beta a partir de R0 estimado
    R0_estimado = 2.5  # Valor inicial (ajustar según necesidad)
    beta = R0_estimado * (gamma + mu)
    
    # Tiempo en días (10 semanas = 70 días)
    t = np.linspace(0, 70, 70)
    
    # Resolver el modelo SIRD
    solucion = odeint(modelo_sird, [S0, I0, R0, D0], t, args=(beta, gamma, mu, N))
    S, I, R, D = solucion.T
    
    # Gráfico
    plt.subplot(2, 2, i+1)
    plt.plot(t, S, colores[i]+'-', label='Susceptibles')
    plt.plot(t, I, colores[i]+'--', label='Infectados')
    plt.plot(t, R, colores[i]+':', label='Recuperados')
    plt.plot(t, D, 'k-', label='Muertos')
    
    # Puntos observados (semanas 1, 2, 3)
    obs_t = [7, 14, 21]
    obs_I = [df.loc[i, 'I0'], 1500 if i==0 else (700 if i==1 else 1200), 
             5000 if i==0 else (2300 if i==1 else 3000)]
    obs_D = [0, 300 if i==0 else (120 if i==1 else 60), 
             650 if i==0 else (1000 if i==1 else 900)]
    
    plt.scatter(obs_t, obs_I, color=colores[i], marker='o', label='Infectados observados')
    plt.scatter(obs_t, obs_D, color='k', marker='x', label='Muertos observados')
    
    plt.title(f"{df.loc[i, 'Ciudad']}\nPoblación: {N:,}\nCFR: {cfr:.1%}")
    plt.xlabel('Días desde inicio')
    plt.ylabel('Número de personas')
    plt.legend()
    plt.grid(True)

plt.tight_layout()
plt.show()

