from models.unifac_method import *
import time
from models.nrtl_method import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from py4j.java_gateway import JavaGateway

# For unifac method, parameters of ddbst were used and they are viewed on: https://www.ddbst.com/published-parameters-unifac.html

# To use the unifac method, it is required to create the components by class Comp(), which must be fed with the groups
# of the desired molecule
# list_group = [(group number, quantity), (group number, quantity)]


def list_fractions(num_points):
    """
    this function creates a list of molar fractions with number of points to the power of two points
    The result list cotains three fractions, the sum of which is equal to 1.
    """
    list_aux = []
    for x1 in [1 - i * (1 / (num_points - 1)) for i in range(num_points)]:
        x2_max = 1 - x1
        for x2 in [i * (x2_max / (num_points - 1)) for i in range(num_points)]:
            x3 = 1 - x1 - x2
            list_aux.append([round(x1, 5), round(x2, 5), round(x3, 5)])
    return list_aux


def NRTL_method(alpha, list_A, list_fraction, Temperature, R):
    dict_aux = {
        'x_1': None,
        'x_2': None,
        'G_mist': None
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

t_0 = time.time()
fractions = list_fractions(num_points=85)
"""
## NRTL Method
list_a = [
    [0, 483.32, 2022],
    [-116.92, 0, 291.60],
    [2382.4, 259.50, 0],
]
alpha_NRTL = 0.2
dict_NRTL = NRTL_method(alpha=alpha_NRTL,
                        list_A=list_a,
                        list_fraction=fractions,
                        Temperature=298.15,
                        R=8.314)



print(t_f-t_0)
ax = plt.axes(projection='3d')

# Data for a three-dimensional line

# Data for three-dimensional scattered points
xdata = dict_NRTL['x_1']
ydata = dict_NRTL['x_2']
zdata = dict_NRTL['G_mist']
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='twilight_shifted')
ax.set_xlabel('$x_1$', fontsize=12)
ax.set_ylabel('$x_2$', fontsize=12)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$G_{mist}/RT$', fontsize=12, rotation=0)
plt.show()

"""
"""
T = (273.15+25)
R = 8.314
x_1 = 0
x_2 = 0
x_3 = 0
comp_1 = Comp(name="acetona", x=x_1, list_group=[(1, 1), (18, 1)])
comp_2 = Comp(name="Acetic acid ethenyl ester", x=x_2, list_group=[(1, 1), (5, 1), (23, 1)])
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
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='twilight_shifted')
ax.set_xlabel('$x_1$', fontsize=12)
ax.set_ylabel('$x_2$', fontsize=12)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$G_{mist}/RT$', fontsize=12, rotation=0)
plt.show()
"""


## COSMO-SAC
list_comps = ['AZOBENZENE', 'ACETONE', 'WATER']

ax = plt.axes(projection='3d')
dict_cosmo = COSMO_SAC_method(list_comps=list_comps,
                              temperature=298.15)
# Data for three-dimensional scattered points
xdata = dict_cosmo['x_1']
ydata = dict_cosmo['x_2']
zdata = dict_cosmo['G_mist']
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='twilight_shifted')
ax.set_xlabel('$x_1$', fontsize=12)
ax.set_ylabel('$x_2$', fontsize=12)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$G_{mist}/RT$', fontsize=12, rotation=0)
plt.show()

t_f = time.time()

print(t_f-t_0)