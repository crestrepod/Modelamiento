#!/usr/bin/env python

# Libraries import
import sys
import os
import json
import csv
import math
import traceback

from matplotlib import pyplot as plt


parameters_file_path    = 'Parameters/parameters.txt'
mesh_file_path          = 'Mesh/'

def save_parameters_file():
    json_body = {
    "model": {
        "T0": T0,
        "T_c": T_c,
        "E_c": E_c,
        "E_e": E_e,
        "x_s": x_s
    },
    "simulation": {
        "t0": t0,
        "t1": t1,
        "n_step": n_step,
        "thickness_cylinder_surface": thickness,
        "thickness_face": thickness_l,
        "export_calculations":export_vals
    },
    "caja_mesh": {
        "width": caja_width,
        "height": caja_height,
        "length": caja_length,
        "mesh_size": caja_msh_s
    },
    "barril_mesh": {
        "diameter": barril_diameter,
        "height": barril_height,
        "mesh_size": barril_msh_s
    },
    "contenedor": contenedor,
    "extracto":extracto,
    "profile_distance":profile_distance,
    "profile_tolerance": profile_tolerance
    }
    with open(parameters_file_path, 'w') as json_file:
        json.dump(json_body,json_file, indent=2)

try:
    # ============================== PROGRAM START ==============================
    print('\nColcafe 2022 - Modelo de enfriamiento de extracto de café')
    print('Presione CTRL+C en cualquier momento para salir')

    # ============================ MODEL PARAMETERS =============================
    if ( os.path.isfile(parameters_file_path) == True ):
        
        with open(parameters_file_path) as json_file:
            parameters = json.load(json_file)
            
            t0          = parameters["simulation"]["t0"]
            t1          = parameters['simulation']['t1']
            n_step      = parameters['simulation']['n_step']
            thickness   = parameters['simulation']['thickness_cylinder_surface']
            thickness_l = parameters['simulation']['thickness_face']
            export_vals = parameters['simulation']['export_calculations']
            
            T0          = parameters['model']['T0']
            T_c         = parameters['model']['T_c']
            E_c           = parameters['model']['E_c']
            E_e       = parameters['model']['E_e']
            x_s        = parameters['model']['x_s']
            
            contenedor  = parameters['contenedor']
            extracto    = parameters['extracto']

            caja_height = parameters['caja_mesh']['height']
            caja_length = parameters['caja_mesh']['length']
            caja_width  = parameters['caja_mesh']['width']
            caja_msh_s  = parameters['caja_mesh']['mesh_size']

            barril_height   = parameters['barril_mesh']['height']
            barril_diameter = parameters['barril_mesh']['diameter']
            barril_msh_s    = parameters['barril_mesh']['mesh_size']
            
            profile_distance    = parameters['profile_distance']
            profile_tolerance   = parameters['profile_tolerance']

    else:
        print('Error: ¡Parámetros del modelo no encontrados!')
        sys.exit(1)
    
    while(1):
        
        print('\n---------- Menu principal ----------')
        print('1. Modificar malla')
        print('2. Ejecutar el modelo')
        print('3. Verificar las regiones del modelo')
        print('4. Generar gráficos')
        print('5. Generar perfiles de temperatura vs posición')
        print('6. Generar perfiles de temperatura vs tiempo')
        print('7. Finalizar programa')
        a = input("Ingrese el número correspondiente:\n")

        if (a == '1'):
            while(1):
                print('---------- Menu de malla -----------')
                print('Contenedor actual:',contenedor)
                if(contenedor == 'barril'):
                    print('Diámetro: \t',barril_diameter)
                    print('Altura: \t',barril_height)
                    
                elif(contenedor == 'caja'):
                    print('Alto: \t',caja_height)
                    print('Largo: \t',caja_length)
                    print('Ancho: \t',caja_width)
                else:
                    print('Error en configuración del contenedor!')
                    sys.exit(1)

                print('1. Cambiar contenedor')
                print('2. Modificar dimensiones')
                print('3. Generar la malla')
                print('4. Regresar al menú principal')
                b = input("Ingrese el número correspondiente:\n")

                if( b == '1'):
                    if (contenedor == 'barril'):
                        contenedor = 'caja'
                    
                    elif(contenedor == 'caja'):
                        contenedor = 'barril'
                    save_parameters_file()

                elif (b == '2'):
                    try:
                        if (contenedor == 'barril'):
                            barril_diameter = float(input("Ingrese el diámetro del barril(m):"))
                            barril_height   = float(input("Ingrese la altura del barril(m):"))
                        
                        elif(contenedor == 'caja'):
                            caja_height = float(input("Ingrese la altura de la caja(m):"))
                            caja_length = float(input("Ingrese la longitud de la caja(m):"))
                            caja_width  = float(input("Ingrese el ancho de la caja(m):"))
                        
                    except KeyboardInterrupt:
                        print('\n¡Ejecución interrumpida!')
                    
                    except:
                        print('Error con las dimensiones ingresadas, intente nuevamente')
                        continue

                    save_parameters_file()
                    print('\nGenerando malla con los nuevos valores...')
                    try:
                        os.system('python geometria.py')
                        os.system('python ./convert_mesh.py -d 3 '+ mesh_file_path+ contenedor +'.mesh ' + mesh_file_path + contenedor +'.vtk')
                        while(1):
                            c = input("¡Proceso exitoso! ¿Visualizar malla generada? (Y/N):\n")

                            if ( c == 'Y' or c == 'y' ):
                                os.system('python ./resview.py '+ mesh_file_path + contenedor +'.vtk -e')
                                break
                            elif (c == 'N' or c=='n' ):
                                break
                            else:
                                print('¡Comando no reconocido!')
                    except:
                        print('¡Error durante la generación de la malla!')
                        sys.exit(1)
                
                elif (b == '3'):
                    print('\nGenerando malla ('+ contenedor +') con las dimensiones actuales...')
                    try:
                        os.system('python geometria.py')
                        os.system('python ./convert_mesh.py -d 3 '+ mesh_file_path + contenedor +'.mesh ' + mesh_file_path + contenedor +'.vtk')
                        while(1):
                            c = input("¡Proceso exitoso! ¿Visualizar malla generada? (Y/N):\n")

                            if ( c == 'Y' or c == 'y' ):
                                os.system('python ./resview.py '+ mesh_file_path + contenedor +'.vtk -e')
                                break
                            elif (c == 'N' or c=='n' ):
                                break
                            else:
                                print('¡Comando no reconocido!')
                    except:
                        print('¡Error durante la generación de la malla!')
                        sys.exit(1)
                    break

                elif (b == '4'):
                    break
                else:
                    print('¡Comando no reconocido!')

        elif (a == '2'):            
            # ======================== RUN SFEPY MODEL ==================================
            print('----------    Modelo    -----------')
            print('Tiempo inicial: ',t0)
            print('Tiempo final: \t',t1)
            print('Número pasos: \t',n_step)
            print('Extracto: \t',extracto)
            a = input("¿Iniciar ejecución del modelo? (Y/N): \n")

            if ( a == 'Y' or a == 'y' ):
                try:
                    os.system('python ./simple.py colcafe.py -o resultados/resultados')
                    print('¡Fin de la ejecución del modelo!')
                    continue
                except:
                    print('¡Error durante la ejecución del modelo!')
                    sys.exit(1)

            elif (a == 'N' or a=='n' ):
                continue
            else:
                print('¡Comando no reconocido!')
                continue
        
        elif (a == '3'):
            while(1):
                c = input("¿Visualizar las regiones del modelo? (Y/N):\n")

                if ( c == 'Y' or c == 'y' ):
                    run_ok = os.system('python ./resview.py ./resultados/resultados_regions.vtk --no-scalar-bars')
                    if( run_ok!=0 ):
                        print('\nError en la visualización. ¿Ya se ejecutó el modelo?')
                    break
                elif (c == 'N' or c=='n' ):
                    break
                else:
                    print('¡Comando no reconocido!')
            continue

        elif (a == '4'):
            print("Generando gráficos a partir de los resultados actuales...")
            print('Tiempo estimado...',str(n_step*6)+'s')
            for a in range(n_step):
                try:
                    print('\nGráfico',str(a+1)+'/'+str(n_step),' - ',str(round( ((a+1)/n_step)*100,2) )+'%')
                    tmp_number = str(a).zfill(len(str(n_step-1)))
                    tmp_file_path = 'resultados/resultados.'+tmp_number+'.vtk'
                    tmp_save_path = 'Graficas_resultados/resultadoT_'+tmp_number+'.png'
                    os.system('python ./postproc.py --3d --scalar-mode cut_plane '+tmp_file_path+' -b --no-show -o '
                            + tmp_save_path + ' --ranges T,'+str(round(T_c))+','+str(round(T0)))
                except:
                    print('¡Error en la generación del gráfico!')
                    sys.exit(1)
            print('\nGraficos generados en /Graficas_resultados')
            continue

        elif (a == '5'):
            print ('Generación de perfiles de temperatura vs posición...')
            print('Tiempo estimado...',str(n_step*15)+'s')
            
            try:

                for a in range(n_step):
                    print('\nPerfil',str(a+1)+'/'+str(n_step),' - ',str(round( ((a+1)/n_step)*100,2) )+'%')
                    tmp_number = str(a).zfill(len(str(n_step-1)))
                    tmp_file_path = 'resultados/resultados.'+tmp_number+'.vtk'
                    tmp_save_path = 'Graficas/'+contenedor+'_'+tmp_number
                    os.system('python ./probe.py colcafe.py '+tmp_file_path+' -o '
                            + tmp_save_path )
                
                print('\nPerfiles generados en /Graficas')
                continue

            except:
                    print('¡Error en la generación de los perfiles!')
                    sys.exit(1)

        elif (a == '6'):
            
            try:    
                print('Generando gráficos de perfiles en el tiempo...')
                profile_out_path_1 = 'Perfiles/perfil_centro.csv'
                profile_out_path_2 = 'Perfiles/perfil_superficie.csv'
                profile_out_path_3 = 'Perfiles/perfil_distancia.csv'

                lines_x = []
                lines_v1 = []
                lines_v2 = []
                lines_v3 = []

                def_pos = 0

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
                                        if math.isclose(float(tmp_s[0]), barril_height/2, rel_tol = profile_tolerance):
                                            lines_v1.append(float(tmp_s[1]))
                                            break
                                    elif( contenedor == 'caja' ):
                                        if math.isclose(float(tmp_s[0]), caja_length/2, rel_tol = profile_tolerance):
                                            lines_v1.append(float(tmp_s[1]))
                                            break

                        # Get superficial profile
                        with open(tmp_file_path) as file_in:
                            
                            pos = 0
                            for line in file_in:
                                tmp_s = line.split(' ')
                                if ( len(tmp_s) == 2 ):
                                    if ( math.isclose(float(tmp_s[1]),T0, rel_tol = profile_tolerance) and def_pos == 0 ):
                                        lines_v2.append(float(tmp_s[1]))
                                        def_pos = pos
                                        break
                                    elif ( def_pos>0 and pos == def_pos ):
                                        lines_v2.append(float(tmp_s[1]))
                                        break
                                pos += 1
                        
                        # Check profile distance does not exceed mesh dimensions
                        if( (contenedor == 'caja' and profile_distance<=caja_length ) or (contenedor == 'barril' and profile_distance<=barril_height ) ):
                            # Get fixed length profile
                            with open(tmp_file_path) as file_in:
                                
                                for line in file_in:
                                    tmp_s = line.split(' ')
                                    if ( len(tmp_s) == 2 ):

                                        if ( math.isclose(float(tmp_s[0]),profile_distance, rel_tol = profile_tolerance )):
                                            lines_v3.append(float(tmp_s[1]))
                                            break
                        else:
                            raise Exception('La distancia para el perfil excede las dimensiones de la malla, verifique el archivo de parámetros')
                
                if( len(lines_x)==0 or len(lines_v2)==0 or len(lines_v3)==0 ):
                    print('No se pudo encontrar los puntos para los perfiles, verifique el valor de tolerancia e intente de nuevo...')
                
                else:

                    print('Generando perfil en punto medio...')
                    with open(profile_out_path_1,'w',newline='') as f1:
                        writer = csv.writer(f1,delimiter =";")
                        writer.writerow(['Step','Tiempo','Temperatura'])
                        for a in range (0, len(lines_x)):
                            row = [a,lines_x[a],lines_v1[a]]
                            writer.writerow(row)
                    # Override decimal delimiter
                    with open(profile_out_path_1,'r') as f1b:
                        data1 = f1b.read()
                        data1 = data1.replace('.', ',')
                    with open(profile_out_path_1,'w') as f1c:
                        f1c.write(data1)

                    plt.figure(1)
                    plt.plot(lines_x,lines_v1)
                    plt.ylim(T_c-10, T0+10)
                    plt.title('Perfil centro de la geometría')
                    plt.xlabel('Tiempo [s]')
                    plt.ylabel('Temperatura [k]')
                    #plt.show()
                    plt.savefig('Perfiles/Perfil_centro.png', bbox_inches='tight')

                    print('Generando perfil en superficie...')
                    with open(profile_out_path_2,'w',newline='') as f2:
                        writer = csv.writer(f2,delimiter =";")
                        writer.writerow(['Step','Tiempo','Temperatura'])
                        for a in range (0, len(lines_x)):
                            row = [a,lines_x[a],lines_v2[a]]
                            writer.writerow(row)
                    # Override decimal delimiter
                    with open(profile_out_path_1,'r') as f2b:
                        data2 = f2b.read()
                        data2 = data2.replace('.', ',')
                    with open(profile_out_path_1,'w') as f2c:
                        f2c.write(data2)

                    plt.figure(2)
                    plt.plot(lines_x,lines_v2)
                    plt.ylim(T_c-10, T0+10)
                    plt.title('Perfil superficie de la geometría')
                    plt.xlabel('Tiempo [s]')
                    plt.ylabel('Temperatura [k]')
                    #plt.show()
                    plt.savefig('Perfiles/Perfil_superficie.png', bbox_inches='tight')

                    print('Generando perfil en distancia fija...')
                    with open(profile_out_path_3,'w',newline='') as f3:
                        writer = csv.writer(f3,delimiter =";")
                        writer.writerow(['Step','Tiempo','Temperatura'])
                        for a in range (0, len(lines_x)):
                            row = [a,lines_x[a],lines_v3[a]]
                            writer.writerow(row)
                    # Override decimal delimiter
                    with open(profile_out_path_1,'r') as f3b:
                        data3 = f3b.read()
                        data3 = data3.replace('.', ',')
                    with open(profile_out_path_1,'w') as f3c:
                        f3c.write(data3)

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
                continue

            except Exception as e:
                    print(traceback.format_exc())
                    print(e)
                    print('¡Error en la generación de los perfiles!')
                    sys.exit(1)

        elif (a == '7'):
            sys.exit(1)
        else:
            print('¡Comando no reconocido!')

except KeyboardInterrupt:
    print('\n¡Ejecución interrumpida!')