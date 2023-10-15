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

        tipo_cava = parameters['Cava']

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
    

        # Propiedades FCS --------------------------------------------------------------------------------------------        
        if extracto == 'STK':
            coeficientes_cp_positivo = [-5.31442424e-08,  3.39414118e-06, -7.63244649e-05,  7.11173325e-04, -1.83933655e-03,  3.34203004e+00]
            coeficientes_cp_negativo = [2.46475116e-05, 3.01248763e-03, 1.39879835e-01, 3.07007061e+00, 3.19581201e+01, 1.32609805e+02]
            
            coeficientes_rho_positivo = [-3.31523246e-03, -5.68521821e-02,  1.12794191e+03]
            coeficientes_rho_negativo = [2.29383704e-04, 2.41273487e-02, 9.16170874e-01, 1.50320862e+01, 1.16052005e+03]
            
            for index1, row in enumerate(T_values):
                for index2, i in enumerate(row):
                    if i > -3.36:
                        cp_array[index1][index2] = sum([coef * i**potencia for potencia, coef in enumerate(reversed(coeficientes_cp_positivo))]) * 1000 *0.8
                        rho_array[index1][index2] = sum([coef * i**potencia for potencia, coef in enumerate(reversed(coeficientes_rho_positivo))]) *1.2
                    else:
                        cp_array[index1][index2] = sum([coef * i**potencia for potencia, coef in enumerate(reversed(coeficientes_cp_negativo))]) * 1000 *1.2
                        rho_array[index1][index2] = sum([coef * i**potencia for potencia, coef in enumerate(reversed(coeficientes_rho_negativo))]) *1.2


             
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
        
        # Propiedades STK --------------------------------------------------------------------------------------------        
        if(extracto=='STK'):      
            coeficientes_k_positivo = [0.0015469,  0.47318634]
            coeficientes_k_negativo = [-2.09751207e-06, -2.23960338e-04, -8.74781185e-03, -1.49616656e-01,  1.10496274e-01]
            
            for index1, row in enumerate(T_values):
                for index2, i in enumerate(row):
                    if i > -3.36:
                        val[index1][index2] = sum([coef * i**potencia for potencia, coef in enumerate(reversed(coeficientes_k_positivo))])
                    else:
                        val[index1][index2] = sum([coef * i**potencia for potencia, coef in enumerate(reversed(coeficientes_k_negativo))]) 
                           


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


      
            #coeficientes función controlador on-off
            a0=-20.31
            a1=-0.4273
            b1=-0.1188
            a2=0.3271
            b2=-0.04853
            a3=-0.2548
            b3=0.1279
            a4=0.1724
            b4=-0.1499
            a5=-0.1052
            b5=0.1352
            a6=-0.04371
            b6=-0.1293
            a7=-0.01711
            b7=0.09
            a8=-0.0001622
            b8=-0.06818
            w1=0.0004193
            
            # for i in range(0, n_step):
            #     step_time = (t1-t0)/(n_step-1)
            #     time_value = i*step_time
            #     print('time_value',time_value)
            time_value=ts.time
            T_cava =273.15+ a0+(a1*nm.cos(time_value*w1))  + (b1*nm.sin(time_value*w1))+ (a2*nm.cos(2*time_value*w1)) +(b2*nm.sin(2*time_value*w1)) +(a3*nm.cos(3*time_value*w1)) +(b3*nm.sin(3*time_value*w1)) +(a4*nm.cos(4*time_value*w1)) +(b4*nm.sin(4*time_value*w1))  +(a5*nm.cos(5*time_value*w1)) +(b5*nm.sin(5*time_value*w1))  +(a6*nm.cos(6*time_value*w1)) +(b6*nm.sin(6*time_value*w1))  +(a7*nm.cos(7*time_value*w1)) +(b7*nm.sin(7*time_value*w1)) +(a8*nm.cos(8*time_value*w1)) +(b8*nm.sin(8*time_value*w1))
        

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
