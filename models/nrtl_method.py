import numpy as np
class nrtl:
    def __init__(self, list_a_ij, fraction, T, R, alpha):
        self.list_a_ij = list_a_ij
        self.fraction = fraction
        self.T = T
        self.R = R
        self.alpha = alpha
        self.tal = []
        self.G = []
        self.phi = []
        self.theta = []
        self.gamma = []

        self.tal_arrang()
        self.g_arrang()
        self.phi_arrang()
        self.theta_arrang()
        self.gamma_arrang()

    def tal_arrang(self):
        for i, I in enumerate(self.list_a_ij):
            self.tal.append([])
            for j, J in enumerate(I):
                self.tal[i].append(J/(self.R*self.T))

    def g_arrang(self):
        for i, tals in enumerate(self.tal):
            self.G.append([])
            for j, tal in enumerate(tals):
                self.G[i].append(np.exp(-self.alpha*tal))

    def phi_arrang(self):
        for i in range(len(self.fraction)):
            phi_aux = 0
            for j in range(len(self.fraction)):
                phi_aux += self.fraction[j]*self.G[j][i]
            self.phi.append(phi_aux)

    def theta_arrang(self):
        for i in range(len(self.fraction)):
            theta_aux = 0
            for j in range(len(self.fraction)):
                theta_aux += self.fraction[j]*self.tal[j][i]*self.G[j][i]
            self.theta.append(theta_aux)

    def gamma_arrang(self):
        for i in range(len(self.fraction)):
            gamma_aux = self.theta[i]/self.phi[i]
            for j in range(len(self.fraction)):
                gamma_aux += self.fraction[j]*self.G[i][j]*(self.tal[i][j]-self.theta[j]/self.phi[j])/self.phi[j]
            self.gamma.append(np.exp(gamma_aux))

