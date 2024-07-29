import os
import openpyxl
from openpyxl.styles import Font, PatternFill
import numpy as np


def encontrar_indice(value, dict):
    for k, v in dict.items():
        if value in v:
            return k
    return None


class SheetsOpener():
    def __init__(self):
        dir_main = (os.getcwd() + '\\models\\sheets\\unifac_main.xlsx')
        dir_sub = (os.getcwd() + '\\models\\sheets\\unifac_sub.xlsx')
        self.main_sheet = openpyxl.load_workbook(dir_main)['main']
        self.sub_sheet = openpyxl.load_workbook(dir_sub)['sub']


# modelo retirado do apendice G van ness

class Comp():
    def __init__(self, name, x, list_group):
        self._sheets_opener = SheetsOpener()
        self.comp = {
            'name': name,
            'x': x,
            'r': 0.0,
            'q': 0.0,
            'J': 0.0,
            'L': 0.0,
            'ln_comb': 0.0,
            'ln_res': 0.0,
            'ln': 0.0,
            'group': None
        }
        self.populate_group(list_group=list_group)

    def __getitem__(self, key):
        return self.comp[key]

    def __setitem__(self, key, value):
        self.comp[key] = value

    def populate_group(self,list_group):
        _dict_aux = {}
        for i in list_group:
            for num_row, row in enumerate(
                    self._sheets_opener.main_sheet.iter_rows(min_row=1,
                                                             max_row=self._sheets_opener.main_sheet.max_row,
                                                             min_col=1,
                                                             max_col=1,
                                                             values_only=True)):
                if i[0] == row[0]:
                    main_group = self._sheets_opener.main_sheet.cell(row=num_row + 1, column=3).value
                    sub_group = self._sheets_opener.main_sheet.cell(row=num_row + 1, column=1).value
                    if not main_group in _dict_aux:
                        _dict_aux[main_group] = { }
                    _R = self._sheets_opener.main_sheet.cell(row=num_row + 1, column=4).value
                    _Q = self._sheets_opener.main_sheet.cell(row=num_row + 1, column=5).value
                    _dict_aux[main_group][sub_group] = {
                        'qtd': i[1],
                        'R': _R,
                        'Q': _Q
                    }
        self.comp['group'] = _dict_aux


class UNIFAC():
    def __init__(self, temp, list_comp):
        self._sheets_opener = SheetsOpener()
        self.comps = {}
        self.temp = temp

        self._soma_R_comps = 0.0
        self._soma_Q_comps = 0.0

        for i, comp in enumerate(list_comp):
            self.comps[f'comp_{i}'] = comp
        self.string_dict = ""

        for i_comp, v_comp in self.comps.items():
            self.string_dict += v_comp['name'] + ';'
        self.string_dict = self.string_dict.rstrip(';')

        self.TE = {self.string_dict: {}
                   }

        self.populate_dict()
        self.ln_gamma_comb()
        self.ln_gamma_res()

    def populate_dict(self):
        for k_comp, v_comp in self.comps.items():
            soma_R = 0.0
            soma_Q = 0.0
            for k_main, v_main in v_comp['group'].items():
                for k_group, v_group in v_main.items():
                    for row in self._sheets_opener.main_sheet.iter_rows(min_row=2,
                                                                        max_row=self._sheets_opener.main_sheet.max_row,
                                                                        min_col=1,
                                                                        max_col=1,
                                                                        values_only=True):
                        if k_group in row:
                            v_group['R'] = self._sheets_opener.main_sheet.cell(row=row[0] + 1, column=4).value
                            v_group['Q'] = self._sheets_opener.main_sheet.cell(row=row[0] + 1, column=5).value
                            soma_R += (v_group['R'] * v_group['qtd'])
                            soma_Q += (v_group['Q'] * v_group['qtd'])
            v_comp['r'] = soma_R
            v_comp['q'] = soma_Q

    def ln_gamma_comb(self):
        self._soma_R_comps = 0.0
        self._soma_Q_comps = 0.0
        for k_comp, v_comp in self.comps.items():
            self._soma_R_comps += v_comp['x'] * v_comp['r']
            self._soma_Q_comps += v_comp['x'] * v_comp['q']

        for k_comp, v_comp in self.comps.items():
            v_comp['J'] = v_comp['r'] / self._soma_R_comps
            v_comp['L'] = v_comp['q'] / self._soma_Q_comps
            _J = v_comp['J']
            _L = v_comp['L']
            _Q = v_comp['q']
            v_comp['ln_comb'] = 1 - _J + np.log(_J) - 5 * _Q * (1 - (_J / _L) + np.log(_J / _L))

    def ln_gamma_res(self):
        unique_groups = set()
        unique_subgroups = set()
        dict_aux = {}
        for k_comp, v_comp in self.comps.items():
            for k_main, v_main in v_comp['group'].items():
                unique_groups.add(k_main)
                if not k_main in dict_aux:
                    dict_aux[k_main] = []
                for k_group in v_main:
                    unique_subgroups.add(k_group)
                    if not k_group in dict_aux[k_main]:
                        dict_aux[k_main].append(k_group)

        orded_subgroup = sorted(list(unique_subgroups))
        orded_group = sorted(list(unique_groups))
        orded_dict = dict(sorted(dict_aux.items()))
        list_aux = []
        for k_i, v_i in orded_dict.items():
            for i in v_i:
                for k_j, v_j in orded_dict.items():
                    for j in v_j:
                        list_aux.append([encontrar_indice(i, orded_dict), encontrar_indice(j, orded_dict)])

        a_mk = np.zeros((len(orded_subgroup), len(orded_subgroup)), dtype=object)
        for i in range(len(list_aux)):
            r = i // len(orded_subgroup)
            c = i % len(orded_subgroup)
            a_mk[r][c] = list_aux[i]

        # gera a matriz de coeficiente combinatorios
        for i in range(len(orded_subgroup)):
            for j in range(len(orded_subgroup)):
                if a_mk[i][j][0] == a_mk[i][j][1]:
                    a_mk[i][j] = 0
                else:
                    if a_mk[i][j][0] < a_mk[i][j][1]:
                        for num_row, row in enumerate(
                                self._sheets_opener.sub_sheet.iter_rows(min_row=self._sheets_opener.sub_sheet.min_row,
                                                                        max_row=self._sheets_opener.sub_sheet.max_row,
                                                                        min_col=1,
                                                                        max_col=2,
                                                                        values_only=True)):
                            if row[0] == a_mk[i][j][0] and row[1] == a_mk[i][j][1]:
                                v = self._sheets_opener.sub_sheet.cell(row=num_row + 1, column=3).value

                        a_mk[i][j] = v
                    else:
                        for num_row, row in enumerate(
                                self._sheets_opener.sub_sheet.iter_rows(min_row=self._sheets_opener.sub_sheet.min_row,
                                                                        max_row=self._sheets_opener.sub_sheet.max_row,
                                                                        min_col=1,
                                                                        max_col=2,
                                                                        values_only=True)):
                            if row[0] == a_mk[i][j][1] and row[1] == a_mk[i][j][0]:
                                v = self._sheets_opener.sub_sheet.cell(row=num_row + 1, column=4).value
                        a_mk[i][j] = v

        e_ki = np.zeros((len(orded_subgroup), len(self.comps)))

        list_aux = np.zeros((len(self.comps), len(orded_subgroup)), dtype=object)

        for num_row, ind in enumerate(self.comps.items()):
            for k_main, v_main in ind[1]['group'].items():
                for k_group, v_group in v_main.items():
                    i = orded_subgroup.index(k_group)
                    list_aux[num_row][i] = v_group['qtd'] * v_group['Q'] / ind[1]['q']

        for num_col, col in enumerate(list_aux):
            for num_row, row in enumerate(col):
                e_ki[num_row][num_col] = row

        # theta
        theta = []
        list_aux = []
        for k_comp, v_comp in self.comps.items():
            list_aux.append(v_comp['x'] * v_comp['q'])
        for i in range(len(unique_subgroups)):
            soma = 0
            for j in range(len(self.comps)):
                soma += e_ki[i][j] * list_aux[j]
            theta.append(soma / self._soma_Q_comps)


        # Temperature dependency
        if self.temp not in self.TE[self.string_dict]:
            self.TE[self.string_dict][self.temp] = []

        tal_mk = np.zeros((len(orded_subgroup), len(orded_subgroup)))
        for num_row, row in enumerate(a_mk):
            for num_col, col in enumerate(row):
                tal_mk[num_row][num_col] = np.exp((-col / self.temp))
        # beta
        beta_ik = np.zeros((len(self.comps), len(orded_subgroup)))

        for i in range(len(self.comps)):
            for k in range(len(unique_subgroups)):
                soma = 0
                for m in range(len(unique_subgroups)):
                    soma += e_ki[m][i] * tal_mk[m][k]
                beta_ik[i][k] = soma

        # sk
        s_k = []
        for i in range(len(unique_subgroups)):
            soma = 0
            for j in range(len(unique_subgroups)):
                soma += theta[j] * tal_mk[j][i]
            s_k.append(soma)

        # calcular o ln de gamma residual
        for i_comp, v_comp in enumerate(self.comps.items()):
            soma = 0
            for k in range(len(unique_subgroups)):
                soma += theta[k] * beta_ik[i_comp][k] / s_k[k] - e_ki[k][i_comp] * np.log(
                    beta_ik[i_comp][k] / s_k[k])
            v_comp[1]['ln_res'] = v_comp[1]['q'] * (1 - soma)
            v_comp[1]['ln'] = np.exp(v_comp[1]['ln_res'] + v_comp[1]['ln_comb'])
            self.TE[self.string_dict][self.temp].append(v_comp[1]['ln'])