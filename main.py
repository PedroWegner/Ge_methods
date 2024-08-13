from models.unifac_method import *
import time
from models.nrtl_method import *
import numpy as np
from py4j.java_gateway import JavaGateway
import matplotlib.pyplot as plt
import sympy as sp
from sympy import hessian, Matrix


# For unifac method, parameters of ddbst were used and they are viewed on: https://www.ddbst.com/published-parameters-unifac.html

# To use the unifac method, it is required to create the components by class Comp(), which must be fed with the groups
# of the desired molecule
# list_group = [(group number, quantity), (group number, quantity)]


def list_fractions(num_points):
    """
    this function creates a list of molar fractions with number of points to the power of two points
    The result list cotains three fractions, the sum of which is equal to 1.
    """
    # list_aux = []
    # for x1 in [1 - i * (1 / (num_points - 1)) for i in range(num_points)]:
    #     x2_max = 1 - x1
    #     for x2 in [i * (x2_max / (num_points - 1)) for i in range(num_points)]:
    #         x3 = 1 - x1 - x2
    #         list_aux.append([round(x1, 5), round(x2, 5), round(x3, 5)])
    list_aux = []
    # Ajustar os incrementos para evitar 0 e 1
    increment = (0.999 - 0.001) / (num_points - 1)

    for x1 in [0.001 + i * increment for i in range(num_points)]:
        x2_max = 1 - x1
        # Ajustar para que x2 também não alcance 0 ou x2_max
        for x2 in [0.001 + i * (x2_max - 0.001) / (num_points - 1) for i in range(num_points)]:
            x3 = 1 - x1 - x2
            # Certificar-se de que x3 também esteja dentro do intervalo desejado
            if 0.001 <= x3 <= 0.999:
                list_aux.append([round(x1, 5), round(x2, 5), round(x3, 5)])
    return list_aux


def NRTL_method(alpha, list_A, list_fraction, Temperature, R):
    dict_aux = {
        'x_1': None,
        'x_2': None,
        'G_mist': None,
        'gammas': None
    }
    gammas_NRTL = []
    g_mist = []
    x_1 = []
    x_2 = []
    for fraction in list_fraction:
        gammas_NRTL.append(nrtl(
            list_a_ij=list_A,
            fraction=fraction,
            T=Temperature,
            R=R,
            alpha=alpha
        ).gamma)

    for i in range(len(fractions)):
        x_1.append(fractions[i][0])
        x_2.append(fractions[i][1])
        g_mist_aux = 0
        for j in range(len(fractions[i])):
            if fractions[i][j] != 0:
                g_mist_aux += fractions[i][j] * np.log(fractions[i][j])  # x_i*ln|x_i|
            g_mist_aux += np.log(gammas_NRTL[i][j]) * fractions[i][j]  # x_i*ln|gamma_i|
        g_mist.append(g_mist_aux)

    dict_aux['x_1'] = x_1
    dict_aux['x_2'] = x_2
    dict_aux['G_mist'] = g_mist
    dict_aux['gammas'] = gammas_NRTL
    return dict_aux


def COSMO_SAC_method(list_comps, temperature):
    gateway = JavaGateway(auto_field=True)
    JCOSMO = gateway.entry_point
    selected_mode = 'COSMO-SAC-HB2 (GAMESS)'
    model = JCOSMO.newModel(selected_mode)
    n_comps = len(list_comps)
    dict_c = {
        'x_1': None,
        'x_2': None,
        'G_mist': None
    }
    comps = gateway.new_array(gateway.jvm.java.lang.String, n_comps)
    x = gateway.new_array(gateway.jvm.double, n_comps)
    for i in range(len(list_comps)):
        comps[i] = list_comps[i]

    # Set components of the mixture
    model.setCompounds(comps)
    # Set temperature of the system
    model.setTemperature(temperature)

    x_1 = []
    x_2 = []
    g_mist = []

    for fraction in fractions:
        x_1.append(fraction[0])
        x_2.append(fraction[1])
        g_mist_aux = 0
        for i, frac in enumerate(fraction):
            if frac != 0:
                g_mist_aux += frac * np.log(frac)
            x[i] = frac
        model.setComposition(x)
        g_mist_aux += model.excessGibbs()
        g_mist.append(g_mist_aux)

    dict_c['x_1'] = x_1
    dict_c['x_2'] = x_2
    dict_c['G_mist'] = g_mist
    return dict_c


"""
T = (273.15+25)
R = 8.314
x_1 = 0
x_2 = 0
x_3 = 0
comp_2 = Comp(name="acetona", x=x_1, list_group=[(1, 1), (18, 1)])
comp_1 = Comp(name="MIC", x=x_2, list_group=[(1, 2), (2, 1), (3, 1), (18, 1)])
comp_3 = Comp(name="water", x=x_3, list_group=[(16, 1)])
eixo_x = []
eixo_y = []
eixo_z = []
for i, fraction in enumerate(fractions):
    x_1 = round(fraction[0], 5)
    x_2 = round(fraction[1], 5)
    x_3 = round(fraction[2], 5)
    comp_1.comp['x'] = x_1
    comp_2.comp['x'] = x_2
    comp_3.comp['x'] = x_3
    gammas = None
    for comps, temperatures in UNIFAC(temp=T, list_comp=[comp_1, comp_2, comp_3]).TE.items():
        for temperature, list_gamma in temperatures.items():
            gammas = np.array(list_gamma, dtype=np.float64)
    ter_1 = 0
    ter_2 = 0
    ter_3 = 0
    print(gammas, x_1, x_2, x_3)
    if x_1 != 0:
        ter_1 = x_1*np.log(x_1)
    if x_2 != 0:
        ter_2 = x_2*np.log(x_2)
    if x_3 != 0:
        ter_3 = x_3*np.log(x_3)
    g_excess = (ter_1 + ter_2 + ter_3)+(x_1*np.log(gammas[0])+x_2*np.log(gammas[1])+x_3*np.log(gammas[2]))
    eixo_x.append(x_1)
    eixo_y.append(x_2)
    eixo_z.append(g_excess)


ax = plt.axes(projection='3d')
# Data for three-dimensional scattered points
xdata = eixo_x
ydata = eixo_y
zdata = eixo_z
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='bone')
ax.set_xlabel('$x_1$', fontsize=12)
ax.set_ylabel('$x_2$', fontsize=12)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$G_{mist}/RT$', fontsize=12, rotation=0)
plt.show()
"""
"""
## COSMO-SAC
list_comps = ['FORMIC_ACID', '1-HEPTANOL', 'WATER']

ax = plt.axes(projection='3d')
dict_cosmo = COSMO_SAC_method(list_comps=list_comps,
                              temperature=700.15)
# Data for three-dimensional scattered points
xdata = dict_cosmo['x_1']
ydata = dict_cosmo['x_2']
zdata = dict_cosmo['G_mist']
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='bone')
ax.set_xlabel('$x_1$', fontsize=12)
ax.set_ylabel('$x_2$', fontsize=12)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$G_{mist}/RT$', fontsize=12, rotation=0)
plt.show()
"""

t_0 = time.time()
# obtainning the fractions list
n = 70
fractions = list_fractions(num_points=n)
## NRTL Method
list_a = [
    [0, -35.285, 258.07],
    [-107.08, 0, -34.93],
    [1342.3, 559.24, 0],
]
temp = 298.65
alpha_NRTL = 0.2
dict_NRTL = NRTL_method(alpha=alpha_NRTL,
                        list_A=list_a,
                        list_fraction=fractions,
                        Temperature=temp,
                        R=1.0)
# determining parameters of the model NRTL
t = {}
for i, items in enumerate(list_a):
    for j in range(len(items)):
        if not f'{i+1}{j+1}' in t:
            t[f'{i+1}{j+1}'] = None
        t[f'{i+1}{j+1}'] = list_a[i][j]/temp

g = {}
for i, items in enumerate(list_a):
    for j in range(len(items)):
        if not f'{i+1}{j+1}' in g:
            g[f'{i+1}{j+1}'] = None
        g[f'{i+1}{j+1}'] = np.exp(-alpha_NRTL*list_a[i][j]/temp)


# assigning G_mixture function
x1, x2 = sp.symbols('x1 x2')
F1 = (x1 * sp.log(x1) + x2 * sp.log(x2) + (1 - x1 - x2) * sp.log(1 - x1 - x2))
F2 = x1*(x1*t['11']*g['11'] + x2*t['21']*g['21'] + (1-x1-x2)*t['31']*g['31'])/(x1*g['11'] + x2*g['21'] + (1-x1-x2)*g['31'])
F3 = x2*(x1*t['12']*g['12'] + x2*t['22']*g['22'] + (1-x1-x2)*t['32']*g['32'])/(x1*g['12'] + x2*g['22'] + (1-x1-x2)*g['32'])
F4 = (1-x1-x2)*(x1*t['13']*g['13'] + x2*t['23']*g['23'] + (1-x1-x2)*t['33']*g['33'])/(x1*g['13'] + x2*g['23'] + (1-x1-x2)*g['33'])
F_G_mix = F1 + F2 + F3 + F4

# derivation G_mixture
dg_dx1 = sp.diff(F_G_mix, x1)
dg_dx2 = sp.diff(F_G_mix, x2)
dg2_dx1x2 = sp.diff(F_G_mix, x1, x2)
dg2_dx1x1 = sp.diff(F_G_mix, x1, x1)
dg2_dx2x2 = sp.diff(F_G_mix, x2, x2)
# obtaining the Hessian matrix of G_mixture
H = hessian(F_G_mix, (x1, x2))
det_H = dg2_dx1x1*dg2_dx2x2-dg2_dx1x2**2
F_K = det_H / (dg_dx1**2 + dg_dx2**2 + 1) ** 2

# starting code lines for plotting
color_line = '#3A3A3A'
# defining functions to plot
G_func = sp.lambdify((x1, x2), F_G_mix, "numpy")
K_func = sp.lambdify((x1, x2), F_K, "numpy")
x1_vals = np.linspace(0.0005, 0.998, 2500)
x2_vals = np.linspace(0.0005, 0.998, 2500)
x1_grid, x2_grid = np.meshgrid(x1_vals, x2_vals)
# creating grids
K_grid = K_func(x1_grid, x2_grid)
G_grid = G_func(x1_grid, x2_grid)
# matrix where K < 0 and K >= 0
G_grid_positive = np.where(K_grid >= 0, G_grid, np.nan)
G_grid_negative = np.where(K_grid < 0, G_grid, np.nan)
# creating the figure
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x1_grid, x2_grid, G_grid_positive, color='slategrey', alpha=0.6) # G where K >= 0
ax.plot_surface(x1_grid, x2_grid, G_grid_negative, color='mediumturquoise', alpha=0.6) # G where K < 0
# plotting the minimized equilibria criterio
ax.plot3D([0.734, 0.006],
          [0.139, 0.019],
          [-0.538934, -0.0690738],
          color=color_line,
          linewidth=2.5,
          linestyle='-',
          alpha=0.9
          )
ax.plot3D([0.589, 0.007],
          [0.247, 0.040],
          [-0.7148538, -0.1126196],
          color=color_line,
          linewidth=2.5,
          linestyle='-',
          alpha=0.9
          )
ax.plot3D([0.513, 0.008],
          [0.299, 0.055],
          [-0.7634141, -0.1356586],
          color=color_line,
          linewidth=2.5,
          linestyle='-',
          alpha=0.9
          )
ax.plot3D([0.427, 0.007],
          [0.353, 0.077],
          [-0.7881635, -0.164091],
          color=color_line,
          linewidth=2.5,
          linestyle='-',
          alpha=0.9
          )
ax.plot3D([0.380, 0.010],
          [0.379, 0.088],
          [-0.7862252, -0.1798107],
          color=color_line,
          linewidth=2.5,
          linestyle='-',
          alpha=0.9
          )

ax.set_xlabel('$x_1$', fontsize=12)
ax.set_ylabel('$x_2$', fontsize=12)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$G_{mist}/RT$', fontsize=12, rotation=0)
# plt.show()
t_f = time.time()
print(t_f - t_0)
