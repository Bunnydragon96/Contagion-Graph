"""
Proyecto SIR Completo con Animación en Vivo por Tipo de Red
Archivo listo para ejecutar: Proyecto_SIR_Completo_Full.py

Características:
- Generación de 3 topologías: Erdős-Rényi, Barabási-Albert, Watts-Strogatz
- Simulación SIR discreta (por ticks) tal como en tu código original
- Cálculo de métricas, centralidades y top-k nodos
- Experimentos con y sin intervenciones (remoción por degree/betweenness)
- Visualizaciones finales (figuras guardadas)
- **Animación en vivo** (una por cada tipo de red) mostrando:
    * la red con nodos coloreados según S/I/R
    * las curvas S(t), I(t), R(t) que se dibujan en tiempo real
- Comentarios explicativos en todo el código

Ejecutar:
    python3 Proyecto_SIR_Completo_Full.py

Nota: Las animaciones abren ventanas matplotlib. Si ejecutas en servidor sin display, usa X forwarding o guarda los frames (no implementado aquí).
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import random
from collections import defaultdict
from matplotlib.animation import FuncAnimation
import os

# -------------------------
# Semilla para reproducibilidad
# -------------------------
RANDOM_SEED = 12345
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# ==========================================================================================
# GENERACIÓN DE GRAFOS (Funciones originales)
# ==========================================================================================

def generar_erdos_renyi(n, p):
    """
    Genera grafo Erdős-Rényi G(n,p)
    Cada par de nodos se conecta con probabilidad p
    """
    G = nx.Graph()
    G.add_nodes_from(range(n))

    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < p:
                G.add_edge(i, j)

    return G


def generar_barabasi_albert(n, m):
    """
    Genera grafo Barabási-Albert
    Cada nuevo nodo se conecta a m nodos existentes con probabilidad
    proporcional a su grado
    """
    G = nx.Graph()

    # Inicializar con m nodos completamente conectados
    m0 = max(m, 2)
    for i in range(m0):
        G.add_node(i)
        for j in range(i):
            G.add_edge(i, j)
    edge_list = []
    for edge in G.edges():
        edge_list.extend(edge)

    for new_node in range(m0, n):
        G.add_node(new_node)
        targets = set()
        attempts = 0
        while len(targets) < m and attempts < 100:
            if edge_list:
                target = random.choice(edge_list)
                if target != new_node:
                    targets.add(target)
            attempts += 1

        for target in targets:
            G.add_edge(new_node, target)
            edge_list.extend([new_node, target])

    return G


def generar_watts_strogatz(n, k, p):
    """
    Genera grafo Watts-Strogatz (small-world)
    1. Crear anillo regular con k vecinos
    2. Recablear cada arista con probabilidad p
    """
    G = nx.Graph()
    G.add_nodes_from(range(n))

    neighbors = k // 2
    for i in range(n):
        for j in range(1, neighbors + 1):
            neighbor = (i + j) % n
            G.add_edge(i, neighbor)

    edges = list(G.edges())
    for u, v in edges:
        if random.random() < p:
            G.remove_edge(u, v)
            attempts = 0
            while attempts < 50:
                new_v = random.randint(0, n - 1)
                if new_v != u and not G.has_edge(u, new_v):
                    G.add_edge(u, new_v)
                    break
                attempts += 1
            if attempts >= 50:
                G.add_edge(u, v)

    return G

# ==========================================================================================
# MÉTRICAS Y UTILIDADES
# ==========================================================================================

def calcular_centralidades(G):
    """Calcula degree centrality y betweenness centrality"""
    degree_cent = nx.degree_centrality(G)
    betweenness_cent = nx.betweenness_centrality(G)

    return degree_cent, betweenness_cent


def obtener_top_k_nodos(G, k, metrica='degree'):
    """
    Obtiene los top-k nodos según la métrica especificada 'degree' o 'betweenness'
    """
    if metrica == 'degree':
        centrality = nx.degree_centrality(G)
    elif metrica == 'betweenness':
        centrality = nx.betweenness_centrality(G)
    else:
        raise ValueError("Métrica debe ser 'degree' o 'betweenness'")

    nodos_ordenados = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
    top_k = [nodo for nodo, _ in nodos_ordenados[:k]]

    return top_k

# ==========================================================================================
# SIMULACIÓN SIR DISCRETA (función original)
# ==========================================================================================

def simular_sir_discreto(G, inicial_infectados, prob_transmision, prob_recuperacion, 
                         max_ticks=200, nodos_removidos=None):
    """
    Simulación SIR discreta por ticks en un grafo G.
    Permite excluir nodos (intervención) mediante 'nodos_removidos' (set).
    """
    n = G.number_of_nodes()
    nodos = list(G.nodes())

    # Excluir nodos removidos (intervención)
    if nodos_removidos:
        nodos = [n for n in nodos if n not in nodos_removidos]

    if len(nodos) == 0:
        return [0], [0], [0], 0, 0

    # Inicializar estados
    estados = {nodo: 'S' for nodo in nodos}

    # Seleccionar infectados iniciales
    infectados_iniciales = random.sample(nodos, min(inicial_infectados, len(nodos)))
    for nodo in infectados_iniciales:
        estados[nodo] = 'I'

    # Historial
    S_hist = []
    I_hist = []
    R_hist = []

    # Simulación por ticks
    for tick in range(max_ticks):
        # Contar estados actuales
        S_count = sum(1 for estado in estados.values() if estado == 'S')
        I_count = sum(1 for estado in estados.values() if estado == 'I')
        R_count = sum(1 for estado in estados.values() if estado == 'R')

        S_hist.append(S_count)
        I_hist.append(I_count)
        R_hist.append(R_count)

        # Condición de término: no hay más infectados
        if I_count == 0:
            break

        # Preparar cambios de estado para este tick
        nuevas_infecciones = []
        nuevas_recuperaciones = []

        # Procesar cada nodo
        for nodo in nodos:
            if estados[nodo] == 'I':
                # Transmisión: intentar infectar vecinos susceptibles
                vecinos = [v for v in G.neighbors(nodo) if v in estados]
                for vecino in vecinos:
                    if estados[vecino] == 'S':
                        if random.random() < prob_transmision:
                            nuevas_infecciones.append(vecino)

                # Recuperacion: con cierta probabilidad
                if random.random() < prob_recuperacion:
                    nuevas_recuperaciones.append(nodo)

        # Aplicar cambios de estado
        for nodo in nuevas_infecciones:
            estados[nodo] = 'I'

        for nodo in nuevas_recuperaciones:
            estados[nodo] = 'R'

    # Calcular métricas
    pico_I = max(I_hist) if I_hist else 0
    tiempo_pico = I_hist.index(pico_I) if I_hist and pico_I > 0 else 0

    return S_hist, I_hist, R_hist, tiempo_pico, pico_I

# ==========================================================================================
# MÉTRICAS DE EVALUACIÓN (función original)
# ==========================================================================================

def calcular_metricas(S_hist, I_hist, R_hist):
    pico_I = max(I_hist) if I_hist else 0
    tiempo_pico = I_hist.index(pico_I) if I_hist and pico_I > 0 else 0
    R_final = R_hist[-1] if R_hist else 0
    duracion = len(I_hist)

    return {
        'pico_I': pico_I,
        'tiempo_pico': tiempo_pico,
        'R_final': R_final,
        'duracion': duracion
    }

# ==========================================================================================
# FUNCIONES NUEVAS: ANIMACIÓN EN TIEMPO REAL DE SIR + RED
# ==========================================================================================

def graficar_sir_y_red_en_vivo(G, historial_estados, title_prefix=""):
    """
    Dibuja una animación donde:
      - La red cambia sus colores según S/I/R.
      - Las curvas S(t), I(t), R(t) se dibujan mientras avanza el tiempo.
    historial_estados: lista de diccionarios {nodo: 'S'/'I'/'R'} por tick
    title_prefix: texto para agregar al título (ej. nombre de la red)
    """
    # Posición fija para consistencia visual
    pos = nx.spring_layout(G, seed=42)

    fig, (ax_red, ax_curvas) = plt.subplots(1, 2, figsize=(14, 6))

    # Listas para las curvas
    S_list, I_list, R_list = [], [], []

    # Inicializar líneas vacías
    linea_S, = ax_curvas.plot([], [], label="Susceptibles (S)", linewidth=2)
    linea_I, = ax_curvas.plot([], [], label="Infectados (I)", linewidth=2)
    linea_R, = ax_curvas.plot([], [], label="Recuperados (R)", linewidth=2)

    ax_curvas.set_xlim(0, max(10, len(historial_estados)))
    ax_curvas.set_ylim(0, max(10, G.number_of_nodes()))
    ax_curvas.set_xlabel("Tiempo (ticks)")
    ax_curvas.set_ylabel("Número de nodos")
    ax_curvas.set_title(f"{title_prefix} - Curvas SIR (en vivo)")
    ax_curvas.legend()

    def actualizar(frame):
        ax_red.clear()
        estados = historial_estados[frame]

        # Conteos por estado
        S = sum(1 for v in estados.values() if v == "S")
        I = sum(1 for v in estados.values() if v == "I")
        R = sum(1 for v in estados.values() if v == "R")

        S_list.append(S)
        I_list.append(I)
        R_list.append(R)

        # Actualizar límites dinámicamente si hace falta
        ax_curvas.set_xlim(0, max(len(S_list), 10))
        ymax = max(G.number_of_nodes(), max(S_list + I_list + R_list) + 1)
        ax_curvas.set_ylim(0, ymax)

        # Actualizar datos de las líneas
        linea_S.set_data(range(len(S_list)), S_list)
        linea_I.set_data(range(len(I_list)), I_list)
        linea_R.set_data(range(len(R_list)), R_list)

        # Asignar colores por estado a los nodos
        colores = []
        for n in G.nodes():
            estado = estados.get(n, "S")
            if estado == "S":
                colores.append("tab:blue")
            elif estado == "I":
                colores.append("tab:red")
            else:
                colores.append("tab:green")

        nx.draw(G, pos, node_color=colores, node_size=80, ax=ax_red, with_labels=False)
        ax_red.set_title(f"{title_prefix} - Red SIR (tick {frame})")

    anim = FuncAnimation(fig, actualizar, frames=len(historial_estados),
                         interval=250, repeat=False)

    plt.tight_layout()
    plt.show()


def simular_sir_en_tiempo_real(G, inicial_infectados=5,
                               prob_transmision=0.3,
                               prob_recuperacion=0.1,
                               max_ticks=150,
                               title_prefix="Red"):
    """
    Simulación SIR que guarda snapshots de estado para animarlos luego.
    Retorna el historial de estados si se necesita.
    """
    nodos = list(G.nodes())
    estados = {n: "S" for n in nodos}

    # Seleccionar infectados iniciales
    infectados_ini = random.sample(nodos, min(inicial_infectados, len(nodos)))
    for n in infectados_ini:
        estados[n] = "I"

    historial = []

    for t in range(max_ticks):
        # Guardar snapshot para animación
        historial.append(estados.copy())

        nuevos_I = []
        nuevos_R = []

        for nodo in nodos:
            if estados[nodo] == "I":
                # Infectar vecinos susceptibles
                for vecino in G.neighbors(nodo):
                    if estados[vecino] == "S" and random.random() < prob_transmision:
                        nuevos_I.append(vecino)
                # Recuperación
                if random.random() < prob_recuperacion:
                    nuevos_R.append(nodo)

        # Aplicar cambios
        for n in nuevos_I:
            estados[n] = "I"
        for n in nuevos_R:
            estados[n] = "R"

        # Terminar si no hay infectados
        if all(estados[n] != "I" for n in nodos):
            historial.append(estados.copy())
            break

    # Ejecutar la función de trazado/animación
    graficar_sir_y_red_en_vivo(G, historial, title_prefix=title_prefix)

    return historial

# ==========================================================================================
# FUNCIONES DE VISUALIZACIÓN FINAL (mantengo las tuyas, mejoradas con comentarios)
# ==========================================================================================

def guardar_figuras_y_mostrar(resultados_promedio):
    """
    Genera y guarda las figuras finales (como en tu código original).
    resultados_promedio: diccionario con las claves por red y escenario.
    """
    nombres_redes = ['Erdős-Rényi', 'Barabási-Albert', 'Watts-Strogatz']
    colores = {'sin_intervencion': 'blue', 'intervencion_degree': 'red', 
               'intervencion_betweenness': 'green'}

    # FIGURA para cada red: S, I, R
    for nombre in nombres_redes:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # S(t)
        ax = axes[0]
        for escenario, color in colores.items():
            res = resultados_promedio[nombre][escenario]
            label = escenario.replace('_', ' ').title()
            ax.plot(res['t'], res['S'], color=color, linewidth=2.5, label=label, alpha=0.8)
        ax.set_xlabel('Tiempo (ticks)')
        ax.set_ylabel('Susceptibles')
        ax.set_title(f'{nombre} - Susceptibles S(t)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # I(t)
        ax = axes[1]
        for escenario, color in colores.items():
            res = resultados_promedio[nombre][escenario]
            ax.plot(res['t'], res['I'], color=color, linewidth=2.5, label=label, alpha=0.8)
        ax.set_xlabel('Tiempo (ticks)')
        ax.set_ylabel('Infectados')
        ax.set_title(f'{nombre} - Infectados I(t)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # R(t)
        ax = axes[2]
        for escenario, color in colores.items():
            res = resultados_promedio[nombre][escenario]
            ax.plot(res['t'], res['R'], color=color, linewidth=2.5, label=label, alpha=0.8)
        ax.set_xlabel('Tiempo (ticks)')
        ax.set_ylabel('Recuperados')
        ax.set_title(f'{nombre} - Recuperados R(t)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        filename = f'fig_{nombre.replace(" ", "_")}_SIR.png'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"✓ Figura guardada: {filename}")
        plt.close(fig)

    # FIGURA: Comparación de picos
    fig4, ax4 = plt.subplots(1, 1, figsize=(10, 6))
    x = np.arange(len(nombres_redes))
    width = 0.25

    picos_sin = [resultados_promedio[n]['sin_intervencion']['pico'] for n in nombres_redes]
    picos_deg = [resultados_promedio[n]['intervencion_degree']['pico'] for n in nombres_redes]
    picos_bet = [resultados_promedio[n]['intervencion_betweenness']['pico'] for n in nombres_redes]

    ax4.bar(x - width, picos_sin, width, label='Sin intervención', alpha=0.7)
    ax4.bar(x, picos_deg, width, label='Intervención (degree)', alpha=0.7)
    ax4.bar(x + width, picos_bet, width, label='Intervención (betweenness)', alpha=0.7)

    ax4.set_xlabel('Red')
    ax4.set_ylabel('Pico de Infectados')
    ax4.set_title('Comparación de Pico de Infectados por Red e Intervención')
    ax4.set_xticks(x)
    ax4.set_xticklabels(nombres_redes)
    ax4.legend()
    plt.tight_layout()
    plt.savefig('fig_comparacion_picos.png', dpi=300, bbox_inches='tight')
    print("✓ Figura guardada: fig_comparacion_picos.png")
    plt.close(fig4)

    # FIGURA: Tamaño final del brote
    fig5, ax5 = plt.subplots(1, 1, figsize=(10, 6))
    R_final_sin = [resultados_promedio[n]['sin_intervencion']['R_final'] for n in nombres_redes]
    R_final_deg = [resultados_promedio[n]['intervencion_degree']['R_final'] for n in nombres_redes]
    R_final_bet = [resultados_promedio[n]['intervencion_betweenness']['R_final'] for n in nombres_redes]

    ax5.bar(x - width, R_final_sin, width, label='Sin intervención', alpha=0.7)
    ax5.bar(x, R_final_deg, width, label='Intervención (degree)', alpha=0.7)
    ax5.bar(x + width, R_final_bet, width, label='Intervención (betweenness)', alpha=0.7)

    ax5.set_xlabel('Red')
    ax5.set_ylabel('Recuperados Finales')
    ax5.set_title('Comparación de Tamaño Final del Brote por Red e Intervención')
    ax5.set_xticks(x)
    ax5.set_xticklabels(nombres_redes)
    ax5.legend()
    plt.tight_layout()
    plt.savefig('fig_comparacion_brote_final.png', dpi=300, bbox_inches='tight')
    print("✓ Figura guardada: fig_comparacion_brote_final.png")
    plt.close(fig5)


# ==========================================================================================
# EXPERIMENTO COMPLETO (INTEGRACIÓN)
# ==========================================================================================

def ejecutar_experimento_completo(n=500, num_simulaciones=10):
    """
    Ejecuta todo el pipeline:
      1) Genera redes
      2) Muestra animación en vivo (una por tipo de red)
      3) Calcula centralidades y top-k
      4) Ejecuta simulaciones con/ sin intervenciones
      5) Promedia resultados y guarda figuras finales
      6) Muestra resumen de métricas
    """
    inicial_infectados = 5
    prob_transmision = 0.3
    prob_recuperacion = 0.1
    k_remover = 20

    print("=" * 80)
    print("SIMULACIÓN SIR EN REDES CON INTERVENCIONES")
    print("=" * 80)
    print(f"Parámetros: n={n}, infectados_ini={inicial_infectados}, beta={prob_transmision}, gamma={prob_recuperacion}")
    print(f"Top-k nodos para remover: {k_remover}")
    print(f"Simulaciones por configuración: {num_simulaciones}")

    # --------------------------
    # PASO 1: Generar redes
    # --------------------------
    print("\n1. GENERANDO REDES")
    p_er = 0.01
    G_ER = generar_erdos_renyi(n, p_er)
    m_ba = 3
    G_BA = generar_barabasi_albert(n, m_ba)
    k_ws = 6
    p_ws = 0.1
    G_WS = generar_watts_strogatz(n, k_ws, p_ws)

    redes = {
        'Erdős-Rényi': G_ER,
        'Barabási-Albert': G_BA,
        'Watts-Strogatz': G_WS
    }

    # Mostrar resumen básico
    for nombre, G in redes.items():
        print(f"\n{nombre}: N={G.number_of_nodes()}, E={G.number_of_edges()}, grado_avg={2*G.number_of_edges()/G.number_of_nodes():.2f}")

    # --------------------------
    # PASO 1.5: Animación en vivo (una por tipo de red)
    # --------------------------
    print("\n" + "="*40)
    print("ANIMACIONES EN TIEMPO REAL (UNA POR TIPO DE RED)")
    print("="*40)
    # Usamos redes más pequeñas para las animaciones para mayor velocidad visual
    small_n = min(200, n)  # límite para animaciones
    # recrear versiones pequeñas de las redes para animar (mismos parámetros)
    G_ER_small = generar_erdos_renyi(small_n, p_er)
    G_BA_small = generar_barabasi_albert(small_n, m_ba)
    G_WS_small = generar_watts_strogatz(small_n, k_ws, p_ws)

    print("\nMostrando animación para Erdős-Rényi (red reducida)...")
    simular_sir_en_tiempo_real(G_ER_small, inicial_infectados=5, prob_transmision=prob_transmision, prob_recuperacion=prob_recuperacion, title_prefix="Erdős-Rényi (pequeña)")

    print("\nMostrando animación para Barabási-Albert (red reducida)...")
    simular_sir_en_tiempo_real(G_BA_small, inicial_infectados=5, prob_transmision=prob_transmision, prob_recuperacion=prob_recuperacion, title_prefix="Barabási-Albert (pequeña)")

    print("\nMostrando animación para Watts-Strogatz (red reducida)...")
    simular_sir_en_tiempo_real(G_WS_small, inicial_infectados=5, prob_transmision=prob_transmision, prob_recuperacion=prob_recuperacion, title_prefix="Watts-Strogatz (pequeña)")

    # --------------------------
    # PASO 2: Calcular centralidades para intervenciones
    # --------------------------
    print("\n2. CALCULANDO CENTRALIDADES")
    top_k_nodos = {}
    for nombre, G in redes.items():
        print(f"\n{nombre}:")
        top_degree = obtener_top_k_nodos(G, k_remover, metrica='degree')
        top_betweenness = obtener_top_k_nodos(G, k_remover, metrica='betweenness')
        top_k_nodos[nombre] = {'degree': set(top_degree), 'betweenness': set(top_betweenness)}
        print(f"  - Top-{k_remover} por degree y betweenness calculados")

    # --------------------------
    # PASO 3: Ejecutar simulaciones (sin intervención y con intervenciones)
    # --------------------------
    print("\n3. EJECUTANDO SIMULACIONES")
    resultados = defaultdict(lambda: defaultdict(list))

    for nombre, G in redes.items():
        print(f"\nProcesando red: {nombre}")
        # Sin intervención
        print("  - Sin intervención")
        for sim in range(num_simulaciones):
            S, I, R, t_pico, pico = simular_sir_discreto(G, inicial_infectados, prob_transmision, prob_recuperacion)
            resultados[nombre]['sin_intervencion'].append({'S': S, 'I': I, 'R': R, 't_pico': t_pico, 'pico': pico})

        # Intervención por degree
        print("  - Intervención: remoción por degree")
        nodos_remover = top_k_nodos[nombre]['degree']
        for sim in range(num_simulaciones):
            S, I, R, t_pico, pico = simular_sir_discreto(G, inicial_infectados, prob_transmision, prob_recuperacion, nodos_removidos=nodos_remover)
            resultados[nombre]['intervencion_degree'].append({'S': S, 'I': I, 'R': R, 't_pico': t_pico, 'pico': pico})

        # Intervención por betweenness
        print("  - Intervención: remoción por betweenness")
        nodos_remover = top_k_nodos[nombre]['betweenness']
        for sim in range(num_simulaciones):
            S, I, R, t_pico, pico = simular_sir_discreto(G, inicial_infectados, prob_transmision, prob_recuperacion, nodos_removidos=nodos_remover)
            resultados[nombre]['intervencion_betweenness'].append({'S': S, 'I': I, 'R': R, 't_pico': t_pico, 'pico': pico})

    # --------------------------
    # PASO 4: Promediar simulaciones
    # --------------------------
    print("\n4. PROMEDIANDO RESULTADOS")
    def promediar_simulaciones(lista_sims):
        max_len = max(len(sim['I']) for sim in lista_sims)
        S_todas, I_todas, R_todas = [], [], []
        for sim in lista_sims:
            S_ext = list(sim['S']) + [sim['S'][-1]] * (max_len - len(sim['S']))
            I_ext = list(sim['I']) + [sim['I'][-1]] * (max_len - len(sim['I']))
            R_ext = list(sim['R']) + [sim['R'][-1]] * (max_len - len(sim['R']))
            S_todas.append(S_ext)
            I_todas.append(I_ext)
            R_todas.append(R_ext)
        return {
            'S': np.mean(S_todas, axis=0),
            'I': np.mean(I_todas, axis=0),
            'R': np.mean(R_todas, axis=0),
            't': np.arange(max_len),
            'pico': np.mean([sim['pico'] for sim in lista_sims]),
            't_pico': np.mean([sim['t_pico'] for sim in lista_sims]),
            'R_final': np.mean([sim['R'][-1] for sim in lista_sims])
        }

    resultados_promedio = {}
    for nombre in redes.keys():
        resultados_promedio[nombre] = {
            'sin_intervencion': promediar_simulaciones(resultados[nombre]['sin_intervencion']),
            'intervencion_degree': promediar_simulaciones(resultados[nombre]['intervencion_degree']),
            'intervencion_betweenness': promediar_simulaciones(resultados[nombre]['intervencion_betweenness'])
        }

    # --------------------------
    # PASO 5: Análisis de métricas y salida por consola
    # --------------------------
    print("\n5. ANALISIS DE METRICAS")
    print(f"\n{'Red':<20} {'Escenario':<30} {'Pico I':<10} {'T Pico':<10} {'R Final':<10} {'Reducción %':<10}")
    print("-"*100)
    for nombre in redes.keys():
        sin_int = resultados_promedio[nombre]['sin_intervencion']
        int_deg = resultados_promedio[nombre]['intervencion_degree']
        int_bet = resultados_promedio[nombre]['intervencion_betweenness']

        print(f"{nombre:<20} {'Sin intervención':<30} {sin_int['pico']:<10.1f} {sin_int['t_pico']:<10.1f} {sin_int['R_final']:<10.1f} {'-':<10}")
        red_deg = ((sin_int['R_final'] - int_deg['R_final']) / sin_int['R_final'] * 100) if sin_int['R_final'] > 0 else 0
        print(f"{'':<20} {'Intervención (degree)':<30} {int_deg['pico']:<10.1f} {int_deg['t_pico']:<10.1f} {int_deg['R_final']:<10.1f} {red_deg:<10.1f}")
        red_bet = ((sin_int['R_final'] - int_bet['R_final']) / sin_int['R_final'] * 100) if sin_int['R_final'] > 0 else 0
        print(f"{'':<20} {'Intervención (betweenness)':<30} {int_bet['pico']:<10.1f} {int_bet['t_pico']:<10.1f} {int_bet['R_final']:<10.1f} {red_bet:<10.1f}")
        print("-"*100)

    # --------------------------
    # PASO 6: Visualizaciones finales (guardar figuras)
    # --------------------------
    print("\n6. GENERANDO Y GUARDANDO VISUALIZACIONES FINALES")
    guardar_figuras_y_mostrar(resultados_promedio)

    print("\nCONCLUSIONES GENERALES:")
    print(" - Se generaron 3 topologías y se compararon bajo SIR.")
    print(" - Se aplicaron intervenciones removiendo top-k nodos por degree/betweenness.")
    print(" - Las figuras han sido guardadas en el directorio actual.")
    print("="*80)


# ==========================================================================================
# EJECUTAR TODO DESDE MAIN
# ==========================================================================================

if __name__ == "__main__":
    # Crear carpeta de salida opcional
    out_dir = os.getcwd()  # usar directorio actual
    print(f"Ejecutando experimento completo. Salida en: {out_dir}")

    # Ejecutar experimento con valores por defecto (ajusta n y num_simulaciones si quieres)
    ejecutar_experimento_completo(n=500, num_simulaciones=10)
