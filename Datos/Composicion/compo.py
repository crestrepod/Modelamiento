# importar librerías
import numpy as np
import json
import matplotlib.pyplot as plt
import pandas as pd


# Obtener y definir las variables de composición 
filename = 'C:\Users\Camilo\Documents\Modelamiento\Colcafe\Colcafe_Camilo\Datos\Composicion\compo.txt'
with open(filename, 'r') as file:
    params = json.load(file)

xs            = params['xs']
humedad       = params['humedad']   
fibra         = params['fibra']   
lipidos       = params['lipidos']   
cenizas       = params['cenizas']   
proteinas     = params['proteinas']   
carbohidratos = params['carbohidratos']  


#CONSTANTES
R  = 8.314                  # Constante Universal de Gases 
T0 = 273.2                  # Temperatura para modelos
Lo = 333.6                  # Calor latente del agua a T0

T = np.array([-40, -35, -30, -25, -20, -12.8, -10, -7, -5, -4.5, -4, -3, 0, 1, 5, 10, 15, 20, 25])     #Temperaturas [°C] 

T_cong          = -32.523 * xs**2 - 0.6781 * xs                         #Temperatura de congelación [°C]

# Fracción de hielo a cada temperatura
x_hielo = np.zeros_like(T)
for n, i in enumerate(T):
    if i < T_cong:
        x_hielo[n] = humedad * (1 - (T_cong / i))

#Fracción de agua  a cada temperatura
x_agua          = humedad - x_hielo                                        


composicion = np.vstack([proteinas      * np.ones_like(T), 
                         lipidos        * np.ones_like(T), 
                         carbohidratos  * np.ones_like(T), 
                         fibra          * np.ones_like(T), 
                         cenizas        * np.ones_like(T), 
                         x_agua, 
                         x_hielo])


#Calculo de Densidad
rho_protein         = 1.3299e3 - (5.1840e-1 * T)
rho_lipidos         = 9.2559e2 - (4.1757e-1 * T)
rho_carbohidratos   = 1.5991e3 - (3.1046e-1 * T)
rho_fibra           = 1.3115e3 - (3.6589e-1 * T)
rho_cenizas         = 2.4238e3 - (2.8063e-1 * T)
rho_agua            = 9.9718e2 + (3.1439e-3 * T) - (3.7574e-3 * T**2)
rho_hielo           = 9.1689e2 - (1.3071e-1 * T)

rho_i = np.vstack([rho_protein, rho_lipidos, rho_carbohidratos, rho_fibra, rho_cenizas, rho_agua, rho_hielo])

# Densidad del extracto
rho_extracto = 1/ (np.sum(composicion / rho_i , axis=0))       

# Densidad T <= T_cong
rho_menor = rho_extracto[T <= T_cong]                       # Toma los valores de densidad donde la Temperatura era menor o igual que la T de congelación
T_menor   = T[T <= T_cong]                                  # Toma los valores de Temperatura que son menores que T de congelación

## Regresión Lineal 
T_model                     = np.linspace(-40,25,50)
grado_rho_menor             = 4                                                 # Grado del polinomio
coeficientes_rho_menor      = np.polyfit(T_menor, rho_menor, grado_rho_menor)   # Realiza una regresión lineal de 4to grado
polinomio_rho_menor         = np.poly1d(coeficientes_rho_menor)
funcion_rho_menor           = polinomio_rho_menor(T_model[T_model <= T_cong])

# Densidad T > T_cong
rho_mayor = rho_extracto[T >= T_cong]                        # Toma los valores de densidad donde la Temperatura era menor o igual que la T de congelación
T_mayor   = T[T >= T_cong]                                   # Toma los valores de Temperatura que son menores que T de congelación

## Regresión Lineal
grado_rho_mayor = 2
coeficientes_rho_mayor      = np.polyfit(T_mayor, rho_mayor, grado_rho_mayor)                # Realiza una regresión lineal de 4to grado
polinomio_rho_mayor         = np.poly1d(coeficientes_rho_mayor)
funcion_rho_mayor           = polinomio_rho_mayor(T_model[T_model >= T_cong])

funcion_rho = np.concatenate([funcion_rho_menor, funcion_rho_mayor])

plt.figure(1)
plt.title('Densidad vs Temperatura')
plt.plot(T, rho_extracto, 'o', label='Densidad calculada')
plt.plot(T_model, funcion_rho, label='Modelo de densidad')
plt.xlabel('Temperatura [°C]')
plt.ylabel(r'Densidad $[kg/m^3]$')
plt.legend()
plt.show(block=False)


xb = 0.4 * proteinas                                                                # Fracción de agua que no se congela
Ms = -(xs * R * T0**2) / ((humedad - xb) * Lo * T_cong)                             # Masa molecular relativa de solidos solubles en la comida

cp_protein         = 2.0082 + (1.2089e-3 * T) - (1.3129e-6 * T**2)
cp_lipidos         = 1.9842 + (1.4733e-3 * T) - (4.8008e-6 * T**2)
cp_carbohidratos   = 1.5488 + (1.9625e-3 * T) - (5.9399e-6 * T**2)
cp_fibra           = 1.8459 + (1.8306e-3 * T) - (4.4609e-6 * T**2)
cp_cenizas         = 1.0926 + (1.8896e-3 * T) - (3.6817e-6 * T**2)
cp_hielo           = 2.0623 + (6.0769e-3 * T) 
cp_agua            = np.where(T <= 0,  4.1289 - (5.3062e-3 * T) + (9.9516e-4 * T**2) , 4.1289 - (9.0864e-5 * T) + (5.4731e-6 * T**2))
cp_matrix = np.vstack([cp_protein, cp_lipidos, cp_carbohidratos, cp_fibra, cp_cenizas, cp_agua, cp_hielo ])

#Capacidad Calorifica del extracto
mask = T < T_cong
cp_extracto = np.zeros_like(T)
cp_extracto[mask] = 1.55 + 1.26 * xs + ( (xs * R * T0**2 ) / (Ms * T[mask]**2 ) )
cp_extracto[~mask] = np.sum(composicion[:, ~mask] * cp_matrix[:, ~mask], axis=0)  

cp_menor = cp_extracto[T < T_cong]    # Toma los valores de densidad donde la Temperatura era menor o igual que la T de congelación

## Regresión Lineal
grado_cp_menor = 5
coeficientes_cp_menor      = np.polyfit(T_menor, cp_menor, grado_cp_menor)                # Realiza una regresión lineal de 4to grado
polinomio_cp_menor         = np.poly1d(coeficientes_cp_menor)
funcion_cp_menor           = polinomio_cp_menor(T_model[T_model < T_cong])

cp_mayor = cp_extracto[T > T_cong]

grado_cp_mayor = 5
coeficientes_cp_mayor      = np.polyfit(T_mayor, cp_mayor, grado_cp_mayor)                # Realiza una regresión lineal de 5to grado
polinomio_cp_mayor         = np.poly1d(coeficientes_cp_mayor)
funcion_cp_mayor           = polinomio_cp_mayor(T_model[T_model > T_cong])

funcion_cp = np.concatenate([funcion_cp_menor, funcion_cp_mayor])

plt.figure(2)
plt.title('Capacidad Calorífica vs Temperatura')
plt.plot(T,cp_extracto, 'o', label='Cp calculada')
plt.plot(T_model, funcion_cp ,label='Modelo de Cp')
plt.xlabel('Temperatura [°C]')
plt.ylabel(r'Cp $[kJ / kg K]$')
plt.legend()
plt.show(block=False)

k_protein         = 1.7881e-1 + (1.1958e-3 * T) - (2.7178e-6 * T**2)
k_lipidos         = 1.8071e-1 - (2.7604e-4 * T) - (1.7749e-7 * T**2)
k_carbohidratos   = 2.0141e-1 + (1.3874e-3 * T) - (4.3312e-6 * T**2)
k_fibra           = 1.8331e-1 + (1.2497e-3 * T) - (3.1683e-6 * T**2)
k_cenizas         = 3.2962e-1 + (1.4011e-3 * T) - (2.9069e-6 * T**2)
k_agua            = 5.7109e-1 + (1.7625e-3 * T) - (6.7036e-6 * T**2)
k_hielo           = 2.2196    - (6.2489e-3 * T) + (1.0154e-4 * T**2)

k_i = np.vstack([k_protein, k_lipidos, k_carbohidratos, k_fibra, k_cenizas, k_agua, k_hielo])
sustancias = {
        "agua"          : {"rho": rho_agua,         "k": k_agua,            "content": x_agua},
        "fibra"         : {"rho": rho_fibra,        "k": k_fibra,           "content": fibra},
        "lipidos"       : {"rho": rho_lipidos,      "k": k_lipidos,         "content": lipidos},
        "cenizas"       : {"rho": rho_cenizas,      "k": k_cenizas,         "content": cenizas},
        "proteinas"     : {"rho": rho_protein,      "k": k_protein,         "content": proteinas},
        "carbohidratos" : {"rho": rho_carbohidratos,"k": k_carbohidratos,   "content": carbohidratos},
        "hielo"         : {"rho": rho_hielo,        "k": k_hielo,           "content": x_hielo}
}

def calculate_xk(sustancias):
    # Definir las sustancias iniciales
    sustancia_base = 'agua'
    rho_agua = sustancias[sustancia_base]['rho']
    k_agua = sustancias[sustancia_base]['k']
    x_agua = sustancias[sustancia_base]['content']
    
    # Diccionario para guardar resultados
    results = {}
    
    # Número de sustancias
    num_sustancias = len(sustancias) - 1
    
    # Inicializar valores
    x_total = x_agua / rho_agua
    k_previous = k_agua
    
    # Iterar sobre las demás sustancias
    for i, sustancia in enumerate(sustancias.keys()):
        # Saltar la sustancia 'agua'
        if sustancia == sustancia_base:
            continue
        
        # Propiedades de la sustancia actual
        rho_current = sustancias[sustancia]['rho']
        k_current = sustancias[sustancia]['k']
        content_current = sustancias[sustancia]['content']
        
        # Calcular x para la sustancia actual
        x_current = content_current / rho_current
        x_total += x_current
        
        # Calcular factor x
        x_factor = x_current / x_total
        
        # Calcular k para la combinación de sustancias
        k_combined = k_previous * (k_current + 2*k_previous - 2*x_factor*(k_previous - k_current)) / \
                     (k_current + 2*k_previous + x_factor*(k_previous - k_current))
        
        # Actualizar k_previous para la siguiente iteración
        k_previous = k_combined
        
        # Guardar resultados
        results[sustancia] = {'x': x_factor, 'k': k_combined}
    
    return results
resultados = calculate_xk(sustancias)

df = pd.DataFrame({(componente, variable): valores for componente, variables in resultados.items() for variable, valores in variables.items()}, index=T)


k_extracto = resultados['hielo']['k']
k_extracto

# Conductividad T <= T_cong
k_menor = k_extracto[T <= T_cong]                       # Toma los valores de densidad donde la Temperatura era menor o igual que la T de congelación
T_menor   = T[T <= T_cong]                                  # Toma los valores de Temperatura que son menores que T de congelación

## Regresión Lineal 
T_model                     = np.linspace(-40,25,50)
grado_k_menor             = 4                                                 # Grado del polinomio
coeficientes_k_menor      = np.polyfit(T_menor, k_menor, grado_k_menor)   # Realiza una regresión lineal de 4to grado
polinomio_k_menor         = np.poly1d(coeficientes_k_menor)
funcion_k_menor           = polinomio_k_menor(T_model[T_model <= T_cong])

# Conductividad T > T_cong
k_mayor = k_extracto[T >= T_cong]                        # Toma los valores de densidad donde la Temperatura era menor o igual que la T de congelación
T_mayor   = T[T >= T_cong]                                   # Toma los valores de Temperatura que son menores que T de congelación

## Regresión Lineal
grado_k_mayor = 1
coeficientes_k_mayor      = np.polyfit(T_mayor, k_mayor, grado_k_mayor)                # Realiza una regresión lineal de 4to grado
polinomio_k_mayor         = np.poly1d(coeficientes_k_mayor)
funcion_k_mayor           = polinomio_k_mayor(T_model[T_model >= T_cong])

funcion_k = np.concatenate([funcion_k_menor, funcion_k_mayor])


plt.figure(3)
plt.title('Conductividad vs Temperatura')
plt.plot(T,k_extracto, 'o', label='k calculada')
plt.plot(T_model, funcion_k ,label='Modelo de k')
plt.xlabel('Temperatura [°C]')
plt.ylabel(r'k')
plt.legend()
plt.show(block=True)


print("Funciones de densidad")
print('T < T_cong: ', coeficientes_rho_menor)
print('T > T_cong: ', coeficientes_rho_mayor)
print('')
print("Funciones de Cp")
print('T < T_cong: ', coeficientes_cp_menor)
print('T > T_cong: ', coeficientes_cp_mayor)
print('')
print("Funciones de k")
print('T < T_cong: ', coeficientes_k_menor)
print('T > T_cong: ', coeficientes_k_mayor)

print('Temperatura Congelación: ', T_cong)