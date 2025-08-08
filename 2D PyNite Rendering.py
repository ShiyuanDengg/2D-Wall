import numpy as np
from Pynite.FEModel3D import FEModel3D
from Pynite import FEModel3D
from Pynite.Material import Material
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

E = 3e10 #Pa
nu = 0.2 # Poisson's ratio
G = E / (2 * (1 + nu)) # Shear modulus
rho = 2400 # kg/m^3

mesh_size = 0.1

width = 5   # m
depth = 8   # m
t = 0.1     # m

x_prop = width / 2

H_exposed = 5.0
D_embed = 3.0
H_total = H_exposed + D_embed
L = H_total
EI = 20e6  # Nm²
gamma = 18.0 * 1e3  # N/m³
Ka = 0.333
Kp_mobilised = 1.5  # Reduced passive pressure coefficient
Kp = 3.0
z_prop = 2.5       # Depth of prop from top [m]
zR_prop = 60e3   # N/m (resisting direction, into wall)
zR_prop *= width

tol = 1e-7

def q_soil_smooth(z, width=0.2):

    # Lines on either side of the break
    m1      =  Ka * gamma                                   # slope for z ≤ H
    m2      = (Ka - Kp_mobilised) * gamma                   # slope for z ≥ H
    c2      =  Kp_mobilised * gamma * H_exposed             # intercept so the
                                                            # two branches meet
    # Smooth switch 0→1 built from tanh
    w = 0.5 * (1 + np.tanh((z - H_exposed) / width))        # 0 left, 1 right

    return (1 - w) * (m1 * z) + w * (m2 * z + c2)

def q_soil(z):
    return np.where(
        z <= H_exposed,
        Ka * gamma * z,
        Ka * gamma * z - Kp_mobilised * gamma * (z - H_exposed)
    )

wall = FEModel3D()

wall.add_material('Concrete', E, G, nu, rho)

#creating meshes
wall.add_rectangle_mesh('MESHES', mesh_size, width, depth, t, 'Concrete', plane = 'XZ', x_control = [width/2], y_control = [z_prop])
wall.meshes['MESHES'].generate()

for node in wall.nodes.values():
  if np.isclose(node.X, width/2, tol) and np.isclose(node.Z, z_prop, tol):
    print(node.X, node.Z)
    wall.add_node_load(node.name, 'FY', zR_prop, 'soil')
    break

for element in list(wall.quads.values()):
    Zavg = (element.i_node.Z + element.j_node.Z + element.m_node.Z + element.n_node.Z)/4
    wall.add_quad_surface_pressure(element.name, q_soil_smooth(Zavg), case='soil')

wall.add_load_combo('soil', {'soil': 1})
wall.merge_duplicate_nodes()

for node in wall.nodes.values():
    if (node.Z <= tol) or (depth - node.Z <= tol):
        wall.def_support(node.name, True, True, True, True, True, True)
    elif (node.X <= tol) or (width - node.X <= tol):
        wall.def_support(node.name, support_RZ = True)

wall.analyze(log=True, check_statics=True)

from Pynite.Rendering import Renderer
renderer = Renderer(wall)
renderer.annotation_size = 0.2
renderer.render_loads = False
renderer.deformed_shape = True
renderer.deformed_scale = 1000
renderer.color_map = 'Qy'
renderer.combo_name = 'soil'
renderer.show_labels = False
renderer.scalar_bar = True
renderer.scalar_bar_text_size = 12
renderer.render_model()

