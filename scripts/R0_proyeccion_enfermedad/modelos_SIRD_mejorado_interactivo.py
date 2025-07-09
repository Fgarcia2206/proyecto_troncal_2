import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from ipywidgets import interact, FloatSlider, Dropdown

# Datos de las ciudades
ciudades = {
    'Norilsk': {'poblacion': 180000, 'infectados_iniciales': 300, 'muertes_semana2': 300, 'muertes_semana3': 650},
    'Kayerkan': {'poblacion': 25000, 'infectados_iniciales': 10, 'muertes_semana2': 120, 'muertes_semana3': 1000},
    'Dudinka': {'poblacion': 20000, 'infectados_iniciales': 20, 'muertes_semana2': 60, 'muertes_semana3': 900}
}

def modelo_sird(y, t, beta, gamma, mu, N):
    S, I, R, D = y
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - (gamma + mu) * I
    dRdt = gamma * I
    dDdt = mu * I
    return [dSdt, dIdt, dRdt, dDdt]

@interact(
    ciudad=Dropdown(options=list(ciudades.keys()), description='Ciudad'),
    R0=FloatSlider(min=1.0, max=5.0, step=0.1, value=2.5, description='R₀'),
    duracion_infeccion=FloatSlider(min=5, max=20, step=1, value=12, description='Duración (días)'),
    CFR=FloatSlider(min=0.1, max=0.8, step=0.05, value=0.45, description='CFR')
)
def simular_epidemia(ciudad, R0, duracion_infeccion, CFR):
    datos = ciudades[ciudad]
    N = datos['poblacion']
    I0 = datos['infectados_iniciales']
    S0 = N - I0
    
    # Cálculo de parámetros
    mu = CFR / duracion_infeccion
    gamma = (1 - CFR) / duracion_infeccion
    beta = R0 * (gamma + mu)
    
    # Simulación
    t = np.linspace(0, 70, 70)
    sol = odeint(modelo_sird, [S0, I0, 0, 0], t, args=(beta, gamma, mu, N))
    S, I, R, D = sol.T
    
    # Puntos observados
    obs_t = [0, 14, 21]  # Días (semana 1, 2, 3)
    obs_I = [I0, datos['infectados_iniciales']*5, datos['infectados_iniciales']*16.7]  # Valores aproximados
    obs_D = [0, datos['muertes_semana2'], datos['muertes_semana3']]
    
    # Visualización
    plt.figure(figsize=(12, 6))
    plt.plot(t, S, 'b-', label='Susceptibles', linewidth=2)
    plt.plot(t, I, 'r--', label='Infectados', linewidth=2)
    plt.plot(t, D, 'k-', label='Muertos', linewidth=2)
    plt.scatter(obs_t, obs_I, color='r', s=100, marker='o', label='Infectados observados')
    plt.scatter(obs_t, obs_D, color='k', s=100, marker='X', label='Muertes observadas')
    
    plt.title(f'Modelo SIRD - {ciudad}\n'
              f'R₀ = {R0:.2f} | CFR = {CFR:.0%} | Duración infección = {duracion_infeccion} días\n'
              f'Población: {N:,} | Proyección muertes semana 10: {int(D[-1]):,} ({D[-1]/N:.1%})',
              pad=20)
    
    plt.xlabel('Días desde el inicio', fontsize=12)
    plt.ylabel('Número de personas', fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    # Resultados numéricos
    print(f"\n🔍 Resultados para {ciudad}:")
    print(f"- Pico de infección: {int(I.max()):,} personas (día {int(t[np.argmax(I)])})")
    print(f"- Total infectados semana 10: {int(I[-1] + R[-1]):,} ({((I[-1]+R[-1])/N):.1%} población)")
    print(f"- Muertes acumuladas semana 10: {int(D[-1]):,} ({D[-1]/N:.1%} población)")