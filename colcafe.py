#!/usr/bin/env python
"""
Heat and Mass transfer modeling for DAN products

"""

from __future__ import absolute_import
from sfepy import data_dir
from sfepy.base.base import output
  
import csv
from pandas import read_csv

import numpy as nm
import json
import os
import sys
import math

# ===============================
#          DEBUG CONTROL     
# ===============================
debug = 0

if(debug):
    import time

# ===============================
#  Variables from parameters file     
# ===============================

file_path = 'Parameters/parameters.txt'

# Check if parameters file exist
if ( os.path.isfile(file_path) == True ):
    
    with open(file_path) as json_file:
        parameters = json.load(json_file)
        
        if(debug):
            print('Reading parameters file:')
            print( parameters )
        
        t0          = parameters["simulation"]["t0"]
        t1          = parameters['simulation']['t1']
        n_step      = parameters['simulation']['n_step']
        thickness   = parameters['simulation']['thickness_cylinder_surface']
        thickness_l = parameters['simulation']['thickness_face']
        export_vals = parameters['simulation']['export_calculations']
        
        T0          = parameters['model']['T0']
        T_c         = parameters['model']['T_c']
        E_c         = parameters['model']['E_c']
        E_e         = parameters['model']['E_e']
        x_s         = parameters['model']['x_s']
        
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
    print('Parameters file not found')
    sys.exit(1)

# ===============================
#      Selección de malla     
# ===============================
if ( contenedor == 'caja' ):
    filename_mesh   = 'Mesh/caja.vtk'
    
elif ( contenedor == 'barril' ):
    filename_mesh   = 'Mesh/barril.vtk'
    barril_radius   = barril_diameter/2

if(debug):
    print("Selected mesh: ", filename_mesh)

# ===============================
#           Funciones            
# ===============================
# ----------------------------------------------------------------------
def get_cylinder_surface(coors, domain=None):
    """
    Selects the remaining surface points around the cylinder
    """
    
    x, y, z = coors[:, 0], coors[:, 1], coors[:, 2]
    limite = barril_radius - thickness
    radio = nm.sqrt(y**2.0 + z**2.0)
    flag = nm.where((radio > limite))[0]
    
    return flag

# ----------------------------------------------------------------------
def get_general_surface(coors, domain=None):
    """
    Selects surface elements for coef3 calculation
    """
    
    x, y, z = coors[:, 0], coors[:, 1], coors[:, 2]
    tf = 2
    btf = 2.5
    if (contenedor == 'barril'):
        limite_inf = barril_radius - btf*thickness
        limite_sup = barril_radius - thickness
        radio = nm.sqrt(y**2.0 + z**2.0)
        flag = nm.where(
                        ( (radio > limite_inf)&(radio < limite_sup) ) &
                        ( (x>thickness_l)                           ) &
                        ( x<(barril_height-thickness_l)             )
                        )[0]
    elif(contenedor == 'caja'):

        flag = nm.where(( 
                        ( (x>thickness_l) & ( x<tf*caja_msh_s )             ) &
                        ( (y>thickness_l) & ( y<(caja_width -thickness_l) ) ) &
                        ( (z>thickness_l) & ( z<(caja_length-thickness_l) ) )
                        )
                        |
                        (
                        ( (x>thickness_l) & ( x<(caja_height -thickness_l) ) ) &
                        ( (y>thickness_l) & ( y<tf*caja_msh_s )              ) &
                        ( (z>thickness_l) & ( z<(caja_length-thickness_l)  ) )
                        )
                        |
                        (
                        ( (x>thickness_l) & ( x<(caja_height -thickness_l) ) ) &
                        ( (y>thickness_l) & ( y<(caja_width -thickness_l)  ) ) &
                        ( (z>thickness_l) & ( z<tf*caja_msh_s )              )
                        )
                        |
                        ( 
                        ( (x<caja_height-thickness_l) & ( x>caja_height-(tf*caja_msh_s) ) ) &
                        ( (y>thickness_l) & ( y<(caja_width -thickness_l)               ) ) &
                        ( (z>thickness_l) & ( z<(caja_length-thickness_l)               ) )
                        )
                        |
                        (
                        ( (x>thickness_l) & ( x<(caja_height -thickness_l)              ) ) &
                        ( (y<caja_width-thickness_l) & ( y>caja_width-(tf*caja_msh_s)   ) ) &
                        ( (z>thickness_l) & ( z<(caja_length-thickness_l)               ) )
                        )
                        |
                        (
                        ( (x>thickness_l) & ( x<(caja_height -thickness_l)              ) ) &
                        ( (y>thickness_l) & ( y<(caja_width -thickness_l)               ) ) &
                        ( (z<caja_length-thickness_l) & ( z>caja_length-(tf*caja_msh_s) ) )
                        )
                        )[0]

        
    return flag


# ----------------------------------------------------------------------
def get_cp_rho(ts, coors, problem, equations=None, mode=None, **kwargs):
    """
    Calculates the coeficient for the 1st heat transfer term
    """
    
    if mode == 'qp':
        
        try:
            T_values = problem.evaluate('ev_volume_integrate.1.Omega(T)', 
                                        mode='qp', verbose=False, copy_materials=True)
            if(debug):
                output('Using T_values')
        except:
            if(debug):
                output('Using T0')
            T_values = nm.ones(coors.shape)*T0
        
        # Convert to °C
        T_values = T_values-273.15
        
        cp_array = nm.ones(coors.shape)
        rho_array = nm.ones(coors.shape)
        
        #print('Type:',type(T_values),"Len: ",len(T_values),'Size',T_values.size)
        #print('Type:',type(coors),"Len: ",len(coors),'Size',coors.size)

        # Propiedades FCS --------------------------------------------------------------------------------------------        
        if (extracto=='FCS'):
            index1 = 0
            for row in T_values:
                index2 = 0
                for value in row:
                    if ( value>-3.6 ):
                        
                                
                        cp_array[index1][index2] =((3.35 + (-1.27E-03*value) + (6.38E-05*(value**2)))*1000)*0.8 #J/kg °C
                                               
                        rho_array[index1][index2] =((-0.112*value) + 1128)*1.2 
                    else:
                       
                                                
                        cp_array[index1][index2] = ((139 + (32.9* value) + (3.1*(value**2)) + (0.139*(value**3)) + (2.97E-03*(value**4)) + (2.41E-05*(value**5)))*1000)*1.2 #J/kg °C
                        
                       
                        
                        rho_array[index1][index2] =((1184 + (20.8*value) + (1.4*(value**2)) + (0.04*(value**3))+ (4.05E-04*(value**4))))*1.2
                    index2 += 1
                
                index1 += 1      
             
        # Propiedades STK --------------------------------------------------------------------------------------------        
        elif (extracto=='STK'):
            #PROPIEDADES A DOS TRAMOS
            index1 = 0
            for row in T_values:
                #print(row)
                index2 = 0
                for value in row:
                    if ( value>-4.8 ):
                        
                        cp00=4.118
                        cp10=-0.003419
                        cp01=-2.56
                        cp20=0.0004177
                        cp11=0.0004968
                        cp02=-0.09021
                        cp30=-9.126e-6
                        cp21=-0.0002003
                        cp12=0.004938  
                        
                        cp_array[index1][index2] =(cp00+(cp10*value)+(cp01*x_s)+(cp20*(value**2))+(cp11*value*x_s)+(cp02*(x_s**2))+(cp30*(value**3))+(cp21*(value**2)*x_s)+(cp12*value*(x_s**2)))*1000
                        #((3.11 + (-1.55E-03*value) + (8.45E-05*(value**2)))*1000)#*0.8 #J/kg °C
                        
                        rho00=1006
                        rho10=-0.08894
                        rho01=342.6
                        rho11=-0.04911
                        rho02=287.1  
                        
                        rho_array[index1][index2] =rho00+(rho10*value)+(rho01*x_s)+(rho11*value*x_s)+(rho02*(x_s**2))
                        #((-0.112*value) + 1176)#*1.2 
                            #otro ajuste de densidad((1176 + (-0.122*value) + (-3.26E-03* (value**2))) )*1.2
                    else:
                        cp00=120.4
                        cp10=-59.26
                        cp01=27.39
                        cp20=147.5
                        cp11=-18.59
                        cp02=2.145
                        cp21=25.23
                        cp12=-0.8877
                        cp03=0.07946
                        cp22=1.013
                        cp13=-0.01143
                        cp04=0.001448
                        cp23=0.01145
                        cp14=-2.642e-6
                        cp05=1.049e-5
                        
                        
                        cp_array[index1][index2] =(cp00+(cp10*x_s)+(cp01*value)+(cp20*(x_s**2))+(cp11*x_s*value)+(cp02*(value**2))+(cp21*(x_s**2)*value)+(cp12*x_s*(value**2))+(cp03*(value**3))+(cp22*(x_s**2)*(value**2))+(cp13*x_s*(value**3))+(cp04*(value**4))+(cp23*(x_s**2)*(value**3))+(cp14*x_s*(value**4))+(cp05*(value**5)))*1000
                        #((114 + (23.8* value) + (2.03*(value**2)) + (0.0845*(value**3)) + (1.69E-03*(value**4)) + (1.3E-05*(value**5)))*1000)#*2.5 #J/kg °C
                        
                        rho00=943.5
                        rho10=5.072
                        rho01=559.3
                        rho20=0.2606
                        rho11=11.39
                        rho02=376.4
                        rho30=0.003623
                        rho21=0.168
                        rho12=-0.06816
                        
                        rho_array[index1][index2] = rho00+(rho10*value)+(rho01*x_s)+(rho20*(value**2))+(rho11*value*x_s)+(rho02*(x_s**2))+(rho30*(value**3))+(rho21*(value**2)*x_s)+(rho12*value*(x_s**2))
                        #((1216 + (10.8*value) + (0.4*(value**2)) + (4.72E-03*(value**3))))#*1.2
                    index2 += 1
                    
                index1 += 1    
            
            
            #PROPIEDADES A 3 TRAMOS
            # index1 = 0
            # for row in T_values:
            #     #print(row)
            #     index2 = 0
            #     for value in row:
            #         if ( value>-4 ):
            #             cp_array[index1][index2] = (((-0.001*value)+3.1221)*1000)*0.8 #J/kg °C
            #             rho_array[index1][index2] = ((-0.7064*value)+1177.7)
            #         elif (value<=-4 and value>-7):
            #             cp_array[index1][index2] =  ((105.06 + (17.812* value) + (0.8335*(value**2)) )*1000)*1.7 #J/kg °C
            #             rho_array[index1][index2] = ((1240.7 + (17.185*value) + (0.7052*(value**2)) ))*1.5
            #         else:
            #             #cp_array[index1][index2] = 190000
            #             cp_array[index1][index2] =  (((0.0756*value) + 5.5227 )*1000)*1.2 #J/kg °C
            #             rho_array[index1][index2] = ((1140.4 + (0.9518*value) + (0.0122*(value**2))))
            #         index2 += 1
                    
            #     index1 += 1   
        
        # Save variables in csv file
        if(export_vals):
            
            file_path = 'Variables/cp_rho_'+str(ts.step) +'.csv'
            print('Exportando variables en',file_path +'...', ' Esto podría tardar...')
            with open(file_path,'w',newline='') as f:
                writer = csv.writer(f,delimiter =";")
                writer.writerow(['X','Y','Z','Temperature','rho','cp'])
                for a in range( len(coors) ):
                    tmp_val = T_values[a][0]
                    row = [coors[a][0],coors[a][1],coors[a][2],tmp_val,rho_array[a][0],cp_array[a][0]]
                    writer.writerow(row)
            # Override decimal delimiter and remove temp []
            with open(file_path,'r') as f:
                data = f.read()
                data = data.replace('.', ',')
                data = data.replace('[', '')
                data = data.replace(']', '')
            with open(file_path,'w') as f:
                f.write(data)

        # Coeficient calculation
        val = cp_array*rho_array
        
        output('rho_cp: min:', val.min(), 'max:', val.max())
        val.shape = (val.shape[0] * val.shape[1], 1, 1)
        
        if(debug):
            output(val)
                
        return {'val' : val}

# ----------------------------------------------------------------------
def get_conductivity(ts, coors, problem, equations=None, mode=None, **kwargs):
    """
    Calculates the coeficient for the 2nd heat transfer term
    """
    if mode == 'qp':
        
        try:
            T_values = problem.evaluate('ev_volume_integrate.1.Omega(T)', 
                                        mode='qp', verbose=False, copy_materials=True)
            if(debug):
                output('Using T_values')
        except:
            if(debug):
                output('Using T0')
            T_values = nm.ones(coors.shape)*T0
            
        # Convert to °C
        T_values = T_values-273.15
        
        val = nm.ones(coors.shape)
        
        # Propiedades FCS --------------------------------------------------------------------------------------------        
        if (extracto=='FCS'):
            index1 = 0
            for row in T_values:
                index2 = 0
                for value in row:
                    if ( value>-3.6 ):
                        val[index1][index2] = ((1.56E-03*value) + 0.473)*0.8
                    else:
                    
                        val[index1][index2] =  ( -8.41E-03 + (-0.172*value) + (-10.4E-03*(value**2)) + (-2.75E-04*(value**3)) + (-2.64E-06*(value**4)))*0.8
                    index2 += 1
                
                index1 += 1
        
        # Propiedades STK --------------------------------------------------------------------------------------------        
        elif(extracto=='STK'):         
            index1 = 0
            for row in T_values:
                index2 = 0
                for value in row:
                    if ( value>-4.8 ):
                        
                        k00=0.5869
                        k10=0.001554
                        k01=-0.3407                       
                    
                        val[index1][index2] =(k00+(k10*value)+(k01*x_s))#*0.8
                        #((1.54E-03*value) + 0.444)
                    else:
                        k00=1.427
                        k10=-4.913
                        k01=-0.1271
                        k20=2.916
                        k11=-0.005876
                        k02=-0.006633
                        k21=0.01018
                        k12=0.0001902
                        k03=-0.0001497
                        k22=0.0002354
                        k13=3.88e-6
                        k04=-1.251e-6
                        
                        
                        val[index1][index2] =(k00+(k10*x_s)+(k01*value)+(k20*(x_s**2))+(k11*x_s*value)+(k02*(value**2))+(k21*(x_s**2)*value)+(k12*x_s*(value**2))+(k03*(value**3))+(k22*(x_s**2)*(value**2))+(k13*x_s*(value**3))+(k04*(value**4)))#*0.8
                        #(-0.0434 + (-0.134*value) + (-7.5E-03*(value**2)) + (-1.85E-04*(value**3)) + (-1.67E-06*(value**4)))
                    index2 += 1
                    
                index1 += 1

        # Save variables in csv file
        if(export_vals):
            
            file_path = 'Variables/conductivity_'+str(ts.step) +'.csv'
            print('Exportando variables en',file_path +'...', ' Esto podría tardar...')
            with open(file_path,'w',newline='') as f:
                writer = csv.writer(f,delimiter =";")
                writer.writerow(['X','Y','Z','Temperature','conductivity'])
                for a in range( len(coors) ):
                    tmp_val = T_values[a][0]
                    row = [coors[a][0],coors[a][1],coors[a][2],tmp_val,val[a][0]]
                    writer.writerow(row)
            # Override decimal delimiter and remove temp []
            with open(file_path,'r') as f:
                data = f.read()
                data = data.replace('.', ',')
                data = data.replace('[', '')
                data = data.replace(']', '')
            with open(file_path,'w') as f:
                f.write(data)
          
        
        output('Conductivity: min:', val.min(), 'max:', val.max())
        val.shape = (val.shape[0] * val.shape[1], 1, 1)
        
        return {'val' : val}

# ----------------------------------------------------------------------
def get_T_c(ts, coors, problem, equations=None, mode=None, **kwargs):
    """
    Calculates the coeficient for the 3rd heat transfer term
    """
    
    if mode == 'qp':
        
        try:
            T_values = problem.evaluate('ev_volume_integrate.1.Omega_surf(T)', 
                                        mode='qp', verbose=False, copy_materials=True)
        
            # T_surf = T0
            # T_cava = T_c 
            if(debug):
                output('Using T_values')
                
            T_surf = T_values.mean()            
                #T_cava = nm.ones(ts.shape)
               
            time_value=ts.time
            descenso = 7000
            constante = 21500
            ciclo = descenso + constante

            
            
            tiempo_mod = time_value % ciclo
            # CICLO DE TEMPERATURA DE LA CAVA
            if tiempo_mod <= descenso:
                pendiente = (26.5 - 19) / descenso              # Pendiente del tramo descendente  
                T_cava = 273.15 + -19 - pendiente * tiempo_mod    
            else:
                T_cava = 273.15 + -26.5   

        except:
            if(debug):
                output('Using T0')
            T_surf = T0
            T_cava = T_c       
             
               
        #CÁLCULO CONVECCIÓN FORZADA 
        T_cf=((T_cava+T_surf)/2)-273.15                     # Temperatura de pelicula °C
        
        #propiedades del aire (Tcava - Text_caja)
        rho_cf=-0.0052*T_cf + 1.2923                        #densidad aire a temperatura de pelicula T_cf [kg/m3]
        mu_cf= ((4.9902E-08)*T_cf) + 1.7293E-05             #viscosidad dinámica aire a temperatura de pelicula T_cf [kg/ms]
        v_cf= mu_cf/rho_cf                                  #viscosidad cinématica aire a temperatura de pelicula T_cf [m2/s]
        cp_cf=(-0.0024*(T_cf**2))- 0.0027*T_cf + 1006.1     # capacidad colorífica aire a temperatura de pelicula T_cf [J/kg °C]
        k_cf=((8E-05)*T_cf) + 0.0236                        # conductividad térmica aire a temperatura de pelicula T_cf [W/m °C]
        
        #Cálculo numero adimensionales
        Re_cf=(barril_diameter*0.04)/v_cf
        Pr_cf=(mu_cf*cp_cf)/k_cf   
        
        #Cálculo coeficiente convección forzada
        Nu_cf=0.683*(Re_cf**0.466)*(Pr_cf**(1/3))
        h_cf=(Nu_cf*k_cf)/barril_diameter
        
        
         #CÁLCULO CONVECCIÓN EN RECINTO CERRADO VERTICAL 
        T_rc=((T_cava+T_surf)/2)-273.15                     # Temperatura de pelicula °C
        
        #propiedades del aire (Tsurf - Tintcaja)
        rho_rc=-0.0052*T_rc + 1.2923                        # densidad aire a temperatura de pelicula T_cf [kg/m3]
        mu_rc= ((4.9902E-08)*T_rc) + 1.7293E-05             # viscosidad dinámica aire a temperatura de pelicula T_cf [kg/ms]
        v_rc= mu_rc/rho_rc                                  # viscosidad cinématica aire a temperatura de pelicula T_cf [m2/s]
        cp_rc=(-0.0024*(T_rc**2))- 0.0027*T_rc + 1006.1     # capacidad colorífica aire a temperatura de pelicula T_cf [J/kg °C]
        k_rc=((8E-05)*T_rc) + 0.0236                        # conductividad térmica aire a temperatura de pelicula T_cf [W/m °C]
        
        #Cálculo numero adimensionales
        Pr_rc=(mu_rc*cp_rc)/k_rc  
        Ra_rc=(9.8*(1/(T_rc+273.15))*abs(T_surf-T_cava)*((0.06)**3)*Pr_rc)/(v_rc**2)
         
        
        #Cálculo coeficiente convección forzada
        Nu_rc=0.18*((Pr_rc/(0.2+Pr_rc))*Ra_rc)**(0.29)
        h_rc=(Nu_rc*k_rc)/barril_diameter
        
        #conductividades térmicas materiales resistencias
        k_e=0.688 #conductividad térmica empaque (vidrio) [W/m °C]
        # k_e=0.34 #conductividad térmica empaque (polietileno de baja densidad) [W/m °C]
        # k_e=45 #conductividad térmica empaque (acero al carbon) [W/m °C]
        k_c=0.05  #conductividad térmica caja (cartón) [W/m °C]
        
        #CÁLCULO COEFICIENTE ENGLOBANTE DE TRANSFERENCIA DE CALOR
       # U=((E_e/k_e)+(1/h_rc)+(E_c/k_c)+(1/h_cf))**(-1)
        U=((E_e/k_e)+(1/h_cf))**(-1)
        #U=0.05
        pre_val = U*(T_surf-T_cava)
        val = nm.ones([1,1,1])*pre_val
        a_time = ts.time
        
        print('preval',pre_val)
        print('U',U) 
        print('Tcava',T_cava)
        print('Tsurf',T_surf) 
        print('ts.time',ts.time)
        print('Pr_rc',Pr_rc)
        print('Ra_rc',Ra_rc)
        print('Nu_rc',Nu_rc)
        print('h_rc',h_rc)
         
        
        if(debug):
            output('Time Step: ', a_time)
            output('Surface temp: ', T_surf)
            output('Cava temp: ', T_cava)
            output('U englobante: ', U)
            output('T_c coeficient: ', pre_val)
        
        return {'val' : val}

# ----------------------------------------------------------------------
def gen_probes(problem):
    # Guía: Examples/linear_elasticity/its2D_4.py
    from sfepy.discrete.probes import LineProbe
    n_point = 1000

    if ( contenedor == 'caja' ):
        p0, p1 = nm.array([caja_height/2, caja_width/2, 0.0]), nm.array([caja_height/2, caja_width/2, caja_length])
    
    elif ( contenedor == 'barril' ):
        p0, p1 = nm.array([0.0, 0.0, 0.0]), nm.array([barril_height, 0.0, 0.0])

    line = LineProbe(p0, p1, n_point, share_geometry=True)
    line.set_options(close_limit=0.5)

    probes = [line]
    labels = ['%s -> %s' % (p0, p1)]

    return probes, labels
    
    
# ----------------------------------------------------------------------
def probe_hook(data, probe, label, problem):
    # Guía: Examples/linear_elasticity/its2D_4.py
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    
    def get_it(name, var_name):
        var = problem.create_variables([var_name])[var_name]
        var.set_data(data[name].data)
        pars, vals = probe(var)
        vals = vals.squeeze()
        return pars, vals
    
    results = {}
    results['T'] = get_it('T', 'T')
    
    fig = plt.figure()
    plt.clf()
    fig.subplots_adjust(hspace=0.4)
    plt.subplot(111)
    pars, vals = results['T']
    plt.plot(pars, vals)    
    plt.ylabel('Temperatura[k]')
    plt.xlabel('Distancia sobre eje X[cm]', fontsize=12)
    plt.title('Perfil de temperatura del modelo', fontsize=14,fontweight='bold')
    plt.ylim([round(T_c-10),round(T0+10)])
    #plt.legend(loc='best', prop=fm.FontProperties(size=10))
    
    return plt.gcf(), results

# ======================================================================
#                       SFEPY PROBLEM DEFINITION                        
# ======================================================================

# ----------------------------------------------------------------------
#(1)    "materials": Propiedades de materiales

material_1 = {
    'name'      : 'coef1',
    'function'  : 'get_cp_rho',
}

material_2 = {
    'name'      : 'coef2',
    'function'  : 'get_conductivity',
}

material_3 = {
    'name'      : 'coef3',
    'function'  : 'get_T_c',
}


# ----------------------------------------------------------------------
#(2)    "regions": Se asignan regiones de la malla

region_01 = {
    'name'      : 'Omega',
    'select'    : 'all',
}

region_02 = {
    'name'      : 'Omega_Left',
    'select'    : 'vertices in (x < %.6f)' %thickness_l,
    'kind'      : 'facet',
}

if (contenedor == 'barril'):
    
    region_03 = {
        'name'      : 'Omega_Right',
        'select'    : 'vertices in (x > %.6f)' %(barril_height-thickness_l),
        'kind'      : 'facet',
    }

    region_04 = {
    'name'      : 'Omega_L',
    'select'    : 'vertices by get_cylinder_surface',
    }

    region_05 = {
    'name'      : 'Omega_surf',
    'select'    : 'vertices by get_general_surface',
}

elif(contenedor == 'caja'):
    
    region_03 = {
    'name'      : 'Omega_Right',
    'select'    : 'vertices in (x > %.6f)' %(caja_height-thickness_l),
    'kind'      : 'facet',
    }
    
    region_04 = {
    'name'      : 'Omega_Top',
    'select'    : 'vertices in (z > %.6f)' %(caja_length-thickness_l),
    'kind'      : 'facet',
    }
    
    region_05 = {
    'name'      : 'Omega_Bottom',
    'select'    : 'vertices in (z < %.6f)' %thickness_l,
    'kind'      : 'facet',
    }
    
    region_06 = {
    'name'      : 'Omega_Back',
    'select'    : 'vertices in (y > %.6f)' %(caja_width-thickness_l),
    'kind'      : 'facet',
    }
    
    region_07 = {
    'name'      : 'Omega_Front',
    'select'    : 'vertices in (y < %.6f)' %thickness_l,
    'kind'      : 'facet',
    }

    region_08 = {
    'name'      : 'Omega_surf',
    'select'    : 'vertices by get_general_surface',
}

# ----------------------------------------------------------------------
#(3)    "fields": definir la aproximación del subdominio

field_1 = {
    'name'          : 'temperature',
    'dtype'         : 'real',
    'shape'         : 'scalar',
    'region'        : 'Omega',
    'approx_order'  : 1,
}


# ----------------------------------------------------------------------
#(4)    "variables": se define la variable desconocida y la de prueba

variable_1 = {
    'name'    : 'T',
    'kind'    : 'unknown field',
    'field'   : 'temperature',
    'order'   : 0,
    'history' : 1,
}

variable_2 = {
    'name'  : 's',
    'kind'  : 'test field',
    'field' : 'temperature',
    'dual'  : 'T',
}


# ----------------------------------------------------------------------
#(5)    "ics":     Condiciones iniciales

# Initial condition for temperature
ic_1 = {
    'name'      : 'ic1',
    'region'    : 'Omega',
    'dofs'      : {'T.0' : T0},
}

# ----------------------------------------------------------------------
#(6)    "ebcs":     Condiciones de frontera

# Boundary conditions for temperature

# ebc_1 = {
#     'name'      : 'ebc1',
#     'region'    : 'Omega_Right',
#     'dofs'      : {'T.0' : T_c},
# }

# ebc_2 = {
#     'name'      : 'ebc2',
#     'region'    : 'Omega_Left',
#     'dofs'      : {'T.0' : T_c},
# }

if (contenedor == 'barril'):
    
  ebc_1 = {
    'name'      : 'ebc1',
    'region'    : 'Omega_Right',
    'dofs'      : {'T.0' : T_c},
}
  

#     ebc_3 = {
#     'name'      : 'ebc3',
#     'region'    : 'Omega_L',
#     'dofs'      : {'T.0' : T_c},
# }

    
elif (contenedor == 'caja'):
    
    # ebc_3 = {
    #     'name'      : 'ebc3',
    #     'region'    : 'Omega_Top',
    #     'dofs'      : {'T.0' : T_c},
    # }

    ebc_4 = {
        'name'      : 'ebc4',
        'region'    : 'Omega_Bottom',
        'dofs'      : {'T.0' : T_c},
    }

    # ebc_5 = {
    #     'name'      : 'ebc5',
    #     'region'    : 'Omega_Front',
    #     'dofs'      : {'T.0' : T_c},
    # }
    
    # ebc_6 = {
    #     'name'      : 'ebc6',
    #     'region'    : 'Omega_Back',
    #     'dofs'      : {'T.0' : T_c},
    # }

# ----------------------------------------------------------------------
#(7.1)    "integrals":  2nd order quadrature over a 3 dimensional space

integral_1 = {
    'name'  : 'i',
    'order' : 1,
}

# ----------------------------------------------------------------------
#(7.2)    "equations":  se define el bloque de ecuaciones.

if (contenedor == 'barril'):

    equations = {
        
        'Heat Transfer' :
        """
        dw_volume_dot.i.Omega( coef1.val, s, dT/dt )
        + dw_laplace.i.Omega( coef2.val, s, T )
        + dw_surface_integrate.i.Omega_Left( coef3.val, s )
        + dw_surface_integrate.i.Omega_Right( coef3.val, s )
        + dw_surface_integrate.i.Omega_L( coef3.val, s )
        = 0  
        """
    }
    
    #   equations = {
        
    #     'Heat Transfer' :
    #     """
    #     dw_volume_dot.i.Omega( coef1.val, s, dT/dt )
    #     + dw_laplace.i.Omega( coef2.val, s, T )
    #     + dw_surface_integrate.i.Omega_Left( coef3.val, s )
    #     = 0  
    #     """
    # }
    
    # equations = {
    #     'Heat Transfer' :
    #     """
    #     dw_volume_dot.i.Omega( coef1.val, s, dT/dt )
    #     + dw_laplace.i.Omega( coef2.val, s, T )
    #     + dw_surface_integrate.i.Omega_L( coef3.val, s )
    #     = 0  
    #     """
    # }
    
elif (contenedor == 'caja'):
    equations = {
        
        'Heat Transfer' :
        """
        dw_volume_dot.i.Omega( coef1.val, s, dT/dt )
        + dw_laplace.i.Omega( coef2.val, s, T )
        + dw_surface_integrate.i.Omega_Left( coef3.val, s )
        + dw_surface_integrate.i.Omega_Right( coef3.val, s )
        + dw_surface_integrate.i.Omega_Top( coef3.val, s )
        + dw_surface_integrate.i.Omega_Bottom( coef3.val, s )
        + dw_surface_integrate.i.Omega_Front( coef3.val, s )
        + dw_surface_integrate.i.Omega_Back( coef3.val, s )
        = 0  
        """
    }
    
# ----------------------------------------------------------------------
#(8.1)    "solvers":    tipo de solver linear y no linear y sus respectivas opciones.
#                       Si no se especifica 'ts', se asume estacionario.

solver_0 = {
    'name'      : 'ls',
    'kind'      : 'ls.scipy_direct',
    #'method'    : 'superlu',
    'use_presolve' : True,
}

solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'         : 20,
    'eps_a'         : 1e-5, 
}

"""
solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    #'i_max'         : 20, #P18
    'i_max'         : 1,
    #'eps_a'         : 1e-13, #P18
    'eps_a'         : 1e-10,
    'eps_r'         : 1.0,
    #'macheps'       : 1e-20, #P18
    'macheps'       : 1e-16,
    #'lin_red'       : 1, # P18
    'lin_red'       : 1e-2, # Linear system error < (eps_a * lin_red). Valor original
    #'lin_red'       : 1e-9, # Linear system error < (eps_a * lin_red).
    #'ls_red'        : 0.01, # P18
    'ls_red'        : 0.1,
    'ls_red_warp'   : 0.001,
    #'ls_on'         : 0.9999, #P18
    'ls_on'         : 1.1,
    #'ls_min'        : 1e-3, #P18
    'ls_min'        : 1e-5,
    'check'         : 0,
    'delta'         : 1e-6,
}
"""


solver_2 = {
    'name' : 'ts',
    'kind' : 'ts.simple',

    't0'            : t0,
    't1'            : t1,
    'dt'            : None,
    'n_step'        : n_step,
    'quasistatic'   : False,
    'verbose'       : 1,
}

# ----------------------------------------------------------------------
#(8.2)    "options":    los solvers a usar se especifican en el bloque de opciones.
#                       Se pueden definir diferentes solvers.

options = {
    'nls'   : 'newton',
    'ls'    : 'ls',
    'ts'    : 'ts',
    'save_times'        : 'all',
    'gen_probes'        : 'gen_probes',
    'probe_hook'        : 'probe_hook',
}

# ----------------------------------------------------------------------
#(9)    "functions": Funciones adicionales implementadas para paŕametros
#                    que varían en el tiempo
functions = {
              'get_cylinder_surface'  :        (get_cylinder_surface,),
              'get_general_surface'   :         (get_general_surface,),
              'get_cp_rho'            :                  (get_cp_rho,),
              'get_conductivity'      :            (get_conductivity,),
              'get_T_c'               :                     (get_T_c,),
              'gen_probes'            :                  (gen_probes,),
              'probe_hook'            :                  (probe_hook,),
                
}
