#!/usr/bin/env python

# Libraries import
import sys
import os
import glob
import subprocess
import json
import csv
import numpy as nm
import math
import platform

from matplotlib import pyplot as plt
from matplotlib.figure import Figure
#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from datetime import datetime
from pandas import read_csv

# ============================ MODEL PARAMETERS =============================
file_path = 'Parameters/parameters.txt'

if ( os.path.isfile(file_path) == True ):
    
    with open(file_path) as json_file:
        parameters = json.load(json_file)
        print('Leyendo parámetros del modelo:')
        
        t0          = parameters["simulation"]["t0"]
        t1          = parameters['simulation']['t1']
        n_step      = parameters['simulation']['n_step']
        
        T0          = parameters['model']['T0']
        T_c         = parameters['model']['T_c']
        W           = parameters['model']['W']
        rho_c       = parameters['model']['rho_c']
        heff        = parameters['model']['heff']
        
        contenedor  = parameters['contenedor']

        print('Tiempo inicial: ',t0)
        print('Tiempo final: \t',t1)
        print('Número pasos: \t',n_step)
        print('Contenedor: \t',contenedor)
else:
    print('¡Parámetros del modelo no encontrados!')

# ======================== GENERATE IMAGES ==================================
n1 = 0
n2 = 0
b = input("¿Generar gráficos a partir de los resultados actuales? (Y/N): \n")

if ( b == 'Y' or b == 'y' ):
    
    n1 = int(input("Step inicial para el procesamiento de los resultados:"))
    n2 = int(input("Step final para el procesamiento de los resultados:"))

    for a in range(n1,n2):
        try:
            tmp_number = str(a).zfill(len(str(n_step-1)))
            tmp_file_path = 'resultados/resultados.'+tmp_number+'.vtk'
            tmp_save_path = 'Graficas_resultados/resultadoT_'+tmp_number+'.png'
            os.system('python ./postproc.py --3d --scalar-mode cut_plane '+tmp_file_path+' -b --no-show -o '
                    + tmp_save_path + ' --ranges T,'+str(round(T_c))+','+str(round(T0)))
        except:
            print('Error en la generación del gráfico')

elif (b == 'N' or b=='n' ):
    print('Continuando sin generar gráficos...')

else:
    print('Comando no reconocido... ¡Ejecución finalizada!')
    sys.exit(1)


# ======================== EXTRACT PROFILES ==================================
print('Generando gráficos de perfiles en el tiempo...')

c = input("¿Generar perfiles de temperatura vs posición a partir de los resultados actuales? (Y/N): \n")

if (c == 'Y' or c == 'y'):

    if (n1 == 0 and n2 == 0):
        n1 = int(input("Step inicial para el procesamiento de los resultados:"))
        n2 = int(input("Step final para el procesamiento de los resultados:"))

    for a in range(n1,n2):
        try:
            tmp_number = str(a).zfill(len(str(n_step-1)))
            tmp_file_path = 'resultados/resultados.'+tmp_number+'.vtk'
            tmp_save_path = 'Graficas/'+contenedor+'_'+tmp_number
            os.system('python ./probe.py colcafe.py '+tmp_file_path+' -o '
                    + tmp_save_path )
        except:
            print('Error en la generación del gráfico')

elif (c == 'N' or c == 'n'):
    print('Continuando sin perfiles vs posición')
else:
    print('Comando no reconocido... ¡Ejecución finalizada!')
    sys.exit(1)

# ======================== SAVE TEMPERATURE DATA ==================================
"""
c = input("¿Generar perfiles de temperatura vs tiempo? (Y/N): \n")

if (c == 'Y' or c == 'y'):

    profile_out_path_1 = 'Perfiles/perfil_centro.csv'
    profile_out_path_2 = 'Perfiles/perfil_superficie.csv'
    profile_out_path_3 = 'Perfiles/perfil_distancia.csv'

    lines_x = []
    lines_v1 = []
    lines_v2 = []
    lines_v3 = []

    caja_height     = 300/1000
    caja_width      = 280/1000
    caja_length     = 310/1000
    barril_height   = 885/1000

    def_pos = 0
    profile_distance = 0.4

    for a in range(0, n_step):
                    
        tmp_file_path = 'graficas/'+ contenedor + '_' + str(a).zfill(len(str(n_step-1))) + '_0.txt'
        # Open results value and get data from center point
        if ( os.path.isfile(tmp_file_path) == True ):
            
            step_time = (t1-t0)/(n_step-1)
            time_value = a*step_time
            lines_x.append( time_value )
            
            # Get middle profile
            with open(tmp_file_path) as file_in:
                
                for line in file_in:
                    tmp_s = line.split(' ')
                    if ( len(tmp_s) == 2 ):

                        if ( contenedor == 'barril' ):
                            if math.isclose(float(tmp_s[0]), barril_height/2, rel_tol = 2e-3):
                                lines_v1.append(float(tmp_s[1]))
                                break
                        elif( contenedor == 'caja' ):
                            if math.isclose(float(tmp_s[0]), math.sqrt(caja_height**2 + caja_width**2 + caja_length**2)/2, rel_tol = 2e-3):
                                lines_v1.append(float(tmp_s[1]))
                                break

            # Get superficial profile
            with open(tmp_file_path) as file_in:
                
                pos = 0
                for line in file_in:
                    tmp_s = line.split(' ')
                    if ( len(tmp_s) == 2 ):

                        if ( math.isclose(float(tmp_s[1]),T0, rel_tol = 2e-3) and def_pos == 0 ):
                            lines_v2.append(float(tmp_s[1]))
                            def_pos = pos 
                            break
                        elif ( def_pos>0 and pos == def_pos ):
                            lines_v2.append(float(tmp_s[1]))
                            break
                    pos += 1
            
            # Get fixed length profile
            with open(tmp_file_path) as file_in:
                
                for line in file_in:
                    tmp_s = line.split(' ')
                    if ( len(tmp_s) == 2 ):

                        if ( math.isclose(float(tmp_s[0]),profile_distance, rel_tol = 2e-3 )):
                            lines_v3.append(float(tmp_s[1]))
                            break


    with open(profile_out_path_1,'w') as f1:
        writer = csv.writer(f1)
        writer.writerow(['Step','Tiempo','Temperatura'])
        for a in range (0, len(lines_x)):
            row = [a,lines_x[a],lines_v1[a]]
            writer.writerow(row)

    plt.figure(1)
    plt.plot(lines_x,lines_v1)
    plt.ylim(T_c-10, T0+10)
    plt.title('Perfil centro de la geometría')
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Temperatura [k]')
    #plt.show()
    plt.savefig('Perfiles/Perfil_centro.png', bbox_inches='tight')

    with open(profile_out_path_2,'w') as f2:
        writer = csv.writer(f2)
        writer.writerow(['Step','Tiempo','Temperatura'])
        for a in range (0, len(lines_x)):
            row = [a,lines_x[a],lines_v2[a]]
            writer.writerow(row)

    plt.figure(2)
    plt.plot(lines_x,lines_v2)
    plt.ylim(T_c-10, T0+10)
    plt.title('Perfil superficie de la geometría')
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Temperatura [k]')
    #plt.show()
    plt.savefig('Perfiles/Perfil_superficie.png', bbox_inches='tight')

    with open(profile_out_path_3,'w') as f3:
        writer = csv.writer(f3)
        writer.writerow(['Step','Tiempo','Temperatura'])
        for a in range (0, len(lines_x)):
            row = [a,lines_x[a],lines_v3[a]]
            writer.writerow(row)

    plt.figure(3)
    plt.plot(lines_x,lines_v3)
    plt.ylim(T_c-10, T0+10)
    plt.title('Perfil distancia fija de la geometría')
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Temperatura [k]')
    #plt.show()
    plt.savefig('Perfiles/Perfil_distancia.png', bbox_inches='tight')

    print('Perfiles generados en:')
    print(profile_out_path_1)
    print(profile_out_path_2)
    print(profile_out_path_3)

else:
    print('¡Ejecución finalizada!')
    sys.exit(1)
"""
print('¡Ejecución finalizada!')