# ==========================
#           pygmsh                        
# ========================== 

import pygmsh
#import meshio
import json
import os

# ============================ MESH PARAMETERS =============================

parameters_file_path = 'Parameters/parameters.txt'
mesh_file_path = 'Mesh/'
if ( os.path.isfile(parameters_file_path) == True ):
        
    with open(parameters_file_path) as json_file:
        parameters = json.load(json_file)
        
        contenedor  = parameters['contenedor']

        caja_height = parameters['caja_mesh']['height']
        caja_length = parameters['caja_mesh']['length']
        caja_width  = parameters['caja_mesh']['width']
        caja_msh_s  = parameters['caja_mesh']['mesh_size']

        barril_height   = parameters['barril_mesh']['height']
        barril_diameter = parameters['barril_mesh']['diameter']
        barril_msh_s    = parameters['barril_mesh']['mesh_size']

if ( contenedor == 'barril' ):
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_max = barril_msh_s
        cyl = geom.add_cylinder([0.0, 0.0, 0.0], [barril_height, 0.0, 0.0], barril_diameter/2)
        mesh = geom.generate_mesh()

elif ( contenedor == 'caja' ):
    with pygmsh.occ.Geometry() as geom:
        geom.characteristic_length_max = caja_msh_s
        cyl = geom.add_box([0.0, 0.0, 0.0], [caja_height, caja_width, caja_length])
        mesh = geom.generate_mesh()

mesh_file_name = mesh_file_path + contenedor + '.mesh'

# Clear old mesh
if ( os.path.isfile(mesh_file_name) == True ):
    os.remove(mesh_file_name)
# Save new mesh
mesh.write(mesh_file_name)