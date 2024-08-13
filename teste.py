import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def F_qeq_expression(ql, Kl):
    qeq = ql*Kl*Ceq / (1 + Kl*Ceq)
    return qeq

def calc_F_qeq(ql, Kl, Ceq_values):
    F_qeq = sp.lambdify(Ceq, F_qeq_expression(ql, Kl), modules='numpy')
    F_qeq_values = F_qeq(Ceq_values)
    return F_qeq_values

# pontos experimentais
qeq_exp = [11.64, 13.23, 15.09, 17.03, 18.07, 18.16]
Ce_exp = [9.87, 16.87, 44.89, 86.75, 195.63, 378.69]

# Defini os simbolos e função
Ceq = sp.symbols('Ceq')
ql = 18.90
Kl = 0.34
# Discretiza
Ceq_values = np.linspace(0, 600, 5000)

# avalia
F_qeq_values = calc_F_qeq(ql=ql, Kl=Kl, Ceq_values=Ceq_values)
F_qeq_values_105 = calc_F_qeq(ql=ql*1.05, Kl=Kl*1.05, Ceq_values=Ceq_values)
F_qeq_values_95 = calc_F_qeq(ql=ql*0.95, Kl=Kl*0.95, Ceq_values=Ceq_values)

#plotagem das curvas
plt.figure(figsize=(8, 6))
plt.plot(Ceq_values, F_qeq_values, label=r'$F_{qeq}(C_{eq})$', color='lime') #curva sem erro
plt.plot(Ceq_values, F_qeq_values_105, label=r'$F_{qeq}(C_{eq})$', color='slategrey') #acrescido 5%
plt.plot(Ceq_values, F_qeq_values_95, label=r'$F_{qeq}(C_{eq})$', color='slategrey') #decrescido 5$
#preenchendo
plt.fill_between(Ceq_values, F_qeq_values_105, F_qeq_values_95, color='mediumturquoise', alpha=0.5)
# plotagem dos pontos
plt.scatter(Ce_exp, qeq_exp, color='black', marker='d', label='experimentais')

# enfeitando
plt.xlabel(r'$C_{eq}$')
plt.ylabel(r'$F_{qeq}$')
plt.title('Para o Matheuszinho')
plt.ylim(0, 25)
plt.grid(False)
plt.show()