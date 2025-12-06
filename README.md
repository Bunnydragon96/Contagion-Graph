# Simulación SIR en Grafos

### Comparación de redes _Erdős-Rényi_, _Barabási-Albert_ y _Watts-Strogatz_

**Autores:** Johnel Cunningham & Luis E. Cardoso  
**Curso:** COMP4017 — Análisis de Algoritmos  
**Universidad de Puerto Rico – Recinto de Mayagüez**  
**Diciembre 2025**

---

# **Dependencias**

Instalar con:

```bash
pip install numpy matplotlib networkx
```
---

# **Cómo Ejecutar el Proyecto**

```bash
git clone https://github.com/Bunnydragon96/Contagion-Graph/

python3 Proyecto_SIR_Completo_Full.py
```

El programa:

1. Genera las 3 redes.

2. Muestra **animaciones SIR en vivo** para cada red.

3. Ejecuta múltiples simulaciones con y sin intervenciones.

4. Calcula **promedios**, picos, tamaños finales del brote, etc.

5. Exporta figuras comparativas.

6. Muestra métricas en consola.


---

# **Parámetros Principales**

Dentro del script:

```python
n = 500                    # Número de nodos
num_simulaciones = 10      # Repeticiones por escenario
inicial_infectados = 5
prob_transmision = 0.3     # β
prob_recuperacion = 0.1    # γ
k_remover = 20             # Nodos eliminados en intervenciones
```

Modificables desde el código en:

```python
ejecutar_experimento_completo(...)
```

---

# **Modelos de Redes Implementados**

|Modelo|Nombre|Parámetros|
|---|---|---|
|ER|Erdős-Rényi|`G(n,p)`|
|BA|Barabási-Albert|`G(n,m)`|
|WS|Watts-Strogatz|`G(n,k,p)`|

Los algoritmos implementan exactamente las definiciones originales:

- **Erdős-Rényi:** aristas independientes con probabilidad _p_
    
- **Barabási-Albert:** apego preferencial en crecimiento
    
- **Watts-Strogatz:** small-world con recableo

---

# **Simulación del Modelo SIR**

Simulación discreta por _ticks_:

- _S_ → Susceptible

- _I_ → Infectado

- _R_ → Recuperado


Incluye:

- contagio por probabilidad β  
- recuperación por probabilidad γ  
- detención automática si I = 0  
- historial S(t), I(t), R(t)

---

# **Estrategias de Intervención**

El sistema implementa:

### Eliminación de nodos por _degree centrality_

### Eliminación por _betweenness centrality_

En ambos casos se remueven los **top-k nodos más importantes**.

---

# **Resultados y Figuras**

El script genera automáticamente:

- S(t), I(t), R(t) por red

- Comparación con y sin intervención

- Comparación de picos

- Tamaño final del brote


Las figuras se guardan como:

```
fig_Erdős-Rényi_SIR.png  
fig_Barabási-Albert_SIR.png  
fig_Watts-Strogatz_SIR.png  
fig_comparacion_picos.png  
fig_comparacion_brote_final.png  
```

---

# **Animación en Vivo**

El proyecto incluye una **animación SIR en tiempo real** para cada topología.

Muestra:

- red con nodos coloreados según S/I/R

- curvas S(t), I(t), R(t) actualizándose dinámicamente

---

# **Escenarios de Simulación**

### **1. Sin intervención**

La red se deja intacta.

### **2. Intervención por Degree**

Se eliminan los nodos con mayor número de conexiones.

### **3. Intervención por Betweenness**

Se eliminan los “puentes críticos” de la red.

### Comparaciones realizadas:

- pico máximo de infectados

- tiempo al pico

- tamaño final del brote

- duración de la epidemia


---

# **Semillas (Reproducibilidad)**

El proyecto usa **una semilla global fija** para generar:

- grafos

- contagios

- recuperaciones

- animaciones


Semilla global definida al inicio:

```python
RANDOM_SEED = 12345
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
```

Esto garantiza:

- resultados reproducibles  
- figuras idénticas entre ejecuciones  

Si quieres cambiar la semilla:

```python
RANDOM_SEED = 777
```
