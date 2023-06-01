import sys
import os
import math
import pandas as pd
import sympy as sp
import numpy as np
from sympy import sympify
import matplotlib.pyplot as plt
from shapely.geometry import LineString

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox

from interface import Ui_MainWindow

pd.options.mode.chained_assignment = None


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4):
        fig = Figure(figsize=(width, height))
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)


class MainWindow:

    # PRINCIPAL
    def __init__(self):

        # INICIALIZAÇÃO
        self.main_win = QMainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.main_win)

        self.ui.stackedWidget.setCurrentWidget(self.ui.home)

        self.ui.mccabe_btn.clicked.connect(lambda: self.ui.stackedWidget.setCurrentWidget(self.ui.mccabe))
        self.ui.home_btn.clicked.connect(lambda: self.ui.stackedWidget.setCurrentWidget(self.ui.home))
        self.ui.fug_btn.clicked.connect(lambda: self.ui.stackedWidget.setCurrentWidget(self.ui.fug))

        # ------------------------ MCCABE INÍCIO -------------------------------
        # Definição padrão dos botões
        self.ui.btn_ideal.setChecked(True)
        self.ui.radio_btn_valor.setChecked(True)
        self.ui.qbox_n_ideal.setVisible(False)
        self.ui.pa.setDisabled(True)
        self.ui.pb.setDisabled(True)

        # Função botões ideal / não ideal
        self.ui.btn_ideal.toggled.connect(lambda: self.ideal_status(self.ui.btn_ideal))
        self.ideal = 'IDEAL'

        # Função botões valor ou pa/pb
        self.ui.radio_btn_valor.toggled.connect(lambda: self.alpha_status(self.ui.radio_btn_valor))
        self.alpha = 'VALOR'

        # Botão de cálculo McCabe clicado
        self.ui.btn_exec.clicked.connect(lambda: self.exec_mccabe(self.ideal))

        # Botão escolher arquivo clicado
        arquivos = list()
        for file in os.listdir('./dados'):
            if file.endswith('.csv'):
                arquivos.append(file.rstrip('.csv'))
        self.ui.qbox_n_ideal.addItems(arquivos)
        self.file_az = f'./dados/{arquivos[0]}.csv'

        self.ui.qbox_n_ideal.activated.connect(self.mist_selec)

        # Salvar imagem gráfico
        self.ui.save_fig.clicked.connect(self.savefig)

        # ------------------------   MCCABE FIM  ----------------------------

        # ------------------------ FUG INÍCIO -------------------------------

        self.ui.comps_num_qBox.addItems(['1', '2'])

        # Quando muda o valor da spinbox, atualiza o valor do combobox pela função
        self.ui.qtd_comps.valueChanged.connect(self.update_table)

        # Lê os dados tabelados e monta uma lista com os nomes
        self.dados_antoine = pd.read_csv('antoine.csv')
        items = list(self.dados_antoine['Substance'])  # Nomes componentes (lista)
        self.ui.comps_qBox.addItems(items)
        self.dados_antoine.set_index('Substance', inplace=True)

        # Seleciona os componentes e guarda nos dicionários
        self.ui.comps_qBox.activated.connect(self.comp_selec)
        self.comps = dict()
        self.key_comps = dict()

        # Determina se o componente é LK ou HK
        self.ui.lk_rdbtn.clicked.connect(self.key_pos)
        self.ui.hk_rdbtn.clicked.connect(self.key_pos)

        # Faz a conexão entre o número e o componente
        self.ui.comps_num_qBox.activated.connect(self.cnum_click)

        # Mudanças na fração molar
        self.ui.mole_frac.valueChanged.connect(self.f_molar)
        self.mole_frac_vls = dict()

        # Faz a pré-visualização dos dados para conferência
        self.ui.preview_btn.clicked.connect(self.check_vls)

        # Executa o cálculo
        self.ui.calc_fug_btn.clicked.connect(self.exec_fug)

        # ------------------------ FUG FIM -------------------------------

    # \/\/ ------------------- FUNÇÕES GERAIS -----------------------\/\/
    def show(self):
        self.main_win.show()

    @staticmethod
    def show_popup(txt):
        msg = QMessageBox()
        msg.setWindowTitle('Erro')
        msg.setText(txt)
        msg.setIcon(QMessageBox.Critical)

        x = msg.exec_()

    @staticmethod
    def check_input(num):
        if ',' in str(num):
            num = num.replace(',', '.')

        try:
            res = float(num)
        except:
            return None

        return res

    # \/\/ ------------------------------- MCCABE INÍCIO ------------------------------- \/\/

    def mist_selec(self):
        name = self.ui.qbox_n_ideal.currentText()
        self.file_az = f'./dados/{name}.csv'

    def ideal_status(self, btn):  # Determina se é ideal ou não
        if btn.isChecked():  # Ideal
            self.ui.qbox_n_ideal.setVisible(False)
            self.ideal = 'IDEAL'

        else:
            self.ui.qbox_n_ideal.setVisible(True)
            self.ideal = 'NAO IDEAL'

    def alpha_status(self, btn):  # Determina se é por valor ou parcial
        if btn.isChecked():  # Valor
            self.ui.pa.setDisabled(True)
            self.ui.pb.setDisabled(True)
            self.ui.alpha_vl.setDisabled(False)
            self.alpha = 'VALOR'
        else:  # Parciais
            self.ui.pa.setDisabled(False)
            self.ui.pb.setDisabled(False)
            self.ui.alpha_vl.setDisabled(True)
            self.alpha = 'PARCIAL'

    def savefig(self):
        try:
            self.canvas.figure.savefig('mccabe_graf.png', dpi=600)
        except:
            return

    def exec_mccabe(self, ideal):

        # Entrada e validação de dados
        xD = round(self.ui.xD.value(), 3)
        xB = round(self.ui.xB.value(), 3)
        xF = round(self.ui.xF.value(), 3)
        n_mur = round(self.ui.n_mur.value(), 2)

        if xD > xF > xB:
            if (0 in [xB, xF, xD]) or (1 in [xB, xF, xD]):
                self.show_popup('Componentes devem ter valor diferente de 0 ou 1!')
                return
            pass
        else:
            self.show_popup('Verificar frações! xD > xF > xB')
            return

        q = self.check_input(self.ui.q_mccabe.text())
        R = self.check_input(self.ui.R.text())
        verificar = [q, R]

        # Se for ideal usa o valor ou calcula a partir de pa/pb
        if ideal == 'IDEAL':
            if self.alpha == 'VALOR':  # Volatilidade relativa (alfa)
                al_value = self.check_input(self.ui.alpha_vl.text())
                verificar.append(al_value)
            else:  # Pa/Pb
                a = self.check_input(self.ui.pa.text())
                b = self.check_input(self.ui.pb.text())
                verificar.extend([a, b])

        else:
            # Se não, utiliza valores presentes no arquivo (não ideal)
            if self.file_az is not None:
                dado_curva = self.file_az
            else:
                return

        # Confere se não há nenhum valor de entrada incorreto ou faltando
        for i in verificar:
            if i is None:
                self.show_popup('Valores incorretos ou faltando')
                return

        # Determina o tipo de dado que vai gerar a curva de equilíbrio
        if ideal == 'IDEAL':
            if self.alpha == 'PARCIAL':
                parcial = (a / b)
                dado_curva = parcial
            elif self.alpha == 'VALOR':
                dado_curva = al_value
            else:
                pass

        self.calc_mccabe(xD, xB, xF, n_mur, q, R, dado_curva)

    def calc_mccabe(self, xD, xB, xF, n_mur, q, R, dado_curva):

        # Canvas no frame do Qt (ajusta e apaga para não ter mais de 1 gráfico por vez)
        if self.ui.horizontalLayout.count() > 2:
            self.ui.horizontalLayout.itemAt(2).widget().setParent(None)

        x_values = np.linspace(0, 1, 51)  # [0, 0.02, 0.04, ...] - Eixo x

        if self.ideal == 'IDEAL':
            alpha = dado_curva
            y_values = (alpha * x_values) / (1 + (alpha - 1) * x_values)  # Valores de y p eq em cond ideais / Eq. (3)

            if n_mur != 1:
                mur_y_values = ((y_values - x_values) * n_mur) + x_values  # Eq. (4)

        else:
            df_mccabe = pd.read_csv(dado_curva)
            try:
                x_values_az = df_mccabe['x1 [mol/mol]'].values
                y_values_az = df_mccabe['y1 [mol/mol]'].values
            except:
                self.show_popup('Formatação arquivo incompatível')
                return

            x0 = x_values_az[0]  # Primeiro ponto: mais a esquerda
            xn = x_values_az[-1]  # Último ponto: mais a direita

            if (xD > xn) or (x0 > xB):
                return

        if q == 1:
            q += 0.0000001  # Eq. (7)

        def intersection(x_vl1, y_vl1, x_vl2, y_vl2):
            line_1 = LineString(np.column_stack((x_vl1, y_vl1)))
            line_2 = LineString(np.column_stack((x_vl2, y_vl2)))

            section = line_1.intersection(line_2)
            x, y = section.xy

            return x[0], y[0]  # (x, y)

        def xy_const(lista, pos):
            # Cria lista de valores constantes: pos = 0 vertical pos = 1 horizontal
            if len(lista) == 1:
                # Se tiver só um item, retorna uma lista com 51 dados de xD (linha reta horizontal) 
                return [xD for x in range(51)]
            else:
                # Se não, retorna uma lista com 51 elems do último item add na lista, para vert/horiz
                return [lista[-1][pos] for x in range(51)]

        def show_stage_num(lista):
            # ESCREVE O NÚMERO DE ESTÁGIOS NO GRÁFICO

            for step, coord in enumerate(lista[1::2]):
                self.canvas.axes.text(coord[0] - 0.025, coord[1] + 0.025, step + 1)

        y_rol = ((R * x_values) / (R + 1)) + (xD / (R + 1))  # VALORES DE Y PARA A RETA RECTIFYING - Eq. (6)
        q_line_y = ((q * x_values) / (q - 1)) - (xF / (q - 1))  # Equação linha q - Eq. (7)

        try:
            coord_intersec = intersection(x_values, y_rol, x_values, q_line_y)
        except IndexError:
            self.show_popup('Além dos limites. Ajustar valores.')
            return

        dados = list()
        dados.append(coord_intersec)

        # Checa se o valor não está acima da curva de equilíbrio
        x1 = [coord_intersec[0], coord_intersec[0]]
        y1 = [0, 1]
        if self.ideal == 'IDEAL':
            y_conf = intersection(x1, y1, x_values, y_values)

            if n_mur != 1:
                y_conf = intersection(x1, y1, x_values, mur_y_values)

        else:
            y_conf = intersection(x1, y1, x_values_az, y_values_az)

        if coord_intersec[1] > y_conf[1]:
            self.show_popup('Valores informados ultrapassam a curva de equilíbrio')
            return
        elif coord_intersec[0] < xB:
            self.show_popup('Valor de q muito baixo')
            return

        x_rect = [dados[0][0], xD]    # valores de x para linha de retificação
        y_rect = [dados[0][1], xD]    # valores de y para linha de retificação
        x_stripp = [xB, dados[0][0]]  # valores de x para linha de esgotamento
        y_stripp = [xB, dados[0][1]]  # valores de y para linha de esgotamento

        # Preparação gráfico
        self.canvas = MplCanvas(self)
        fig = plt.figure(num=None, figsize=(7, 6), dpi=600)

        # Início plotagem + cálculos
        if self.ideal == 'IDEAL':
            self.canvas.axes.plot(x_values, y_values, color='r')  # CURVA DE EQUILÍBRIO

            if n_mur != 1:
                self.canvas.axes.plot(x_values, mur_y_values, color='g', linestyle='--')

        else:
            self.canvas.axes.plot(x_values_az, y_values_az, color='r')  # CURVA DE NÃO IDEAL

        self.canvas.axes.plot(x_values, x_values, color='k')  # CURVA 45°

        # DADOS[0][0] = x || DADOS[0][1] = y || intersecção ROL e q-line
        self.canvas.axes.plot(x_rect, y_rect, color='k', linestyle='--')  # RETIFICAÇÃO
        self.canvas.axes.plot([dados[0][0], xF], [dados[0][1], xF], color='k', linestyle='-.')  # ALIMENTAÇÃO
        self.canvas.axes.plot(x_stripp, y_stripp, color='k', linestyle=':')  # ESGOTAMENTO/STRIPPIN

        count = 1
        stage = 0
        feed_stage = 0

        while True:

            # PAR - x constante / reta vertical
            if count % 2 == 0:
                stage += 1

                # VARIÁVEL CONSTANTE (AQUI X=0) (x,y)
                dado_const = xy_const(dados, 0)

                if (dados[-1][0]) <= xB:  # SE O VALOR FOR MENOR QUE xB GERA A RETA VERTICAL E ENCERRA
                    self.canvas.axes.plot([dados[-1][0], dados[-1][0]], [dados[-1][0], dados[-1][1]], color='b')
                    # plt.plot(X DO ULTIMO VALOR, X DO ULTIMO VALOR], [X DO ULTIMO VALOR, Y DO ULTIMO VALOR])
                    break

                # SE O ÚLTIMO VALOR FOR MENOR DO QUE A INTERSECÇÃO ENTRE Q E RETIFICAÇÃO, DESCE ATÉ ESGOTAMENTO
                elif (dados[-1][0]) <= (dados[0][0]):
                    # ESGOTAMENTO
                    coord_intersec = intersection(dado_const, x_values, x_stripp, y_stripp)

                    if feed_stage == 0:
                        feed_stage = stage

                else:
                    # RETIFICAÇÃO
                    coord_intersec = intersection(dado_const, x_values, x_rect, y_rect)

                dados.append(coord_intersec)
                self.canvas.axes.plot([dados[-2][0], dados[-1][0]], [dados[-2][1], dados[-1][1]], color='b')
                # PLOTA VALORES DE X E Y VINDOS DE CIMA ATÉ A INTERSECÇÃO

                count += 1

            # ÍMPAR - y constante / reta horizontal
            else:

                # VARIÁVEL CONSTANTE (AQUI Y=1) (x,y)
                dado_const = xy_const(dados, 1)

                x_values_copy = x_values.copy()
                if self.ideal != 'IDEAL':
                    x_values_copy = x_values_az
                    y_values = y_values_az
                elif n_mur != 1:
                    y_values = mur_y_values

                # COMEÇA AQUI !!
                if count == 1:

                    # INTERSECÇÃO ENTRE CURVA DE EQ E ESTÁGIO
                    coord_intersec = intersection(x_values, dado_const, x_values_copy, y_values)
                    dados.append(coord_intersec)

                    # VAI DO XD, XD ATÉ A CURVA EQUILÍBRIO
                    self.canvas.axes.plot([dados[-1][0], xD], [dados[-1][1], xD], color='b')

                else:
                    coord_intersec = intersection(x_values, dado_const, x_values_copy, y_values)
                    dados.append(coord_intersec)
                    self.canvas.axes.plot([dados[-2][0], dados[-1][0]], [dados[-2][1], dados[-1][1]], color='b')
                    # PLOTA VALORES DE X E Y VINDOS DO LADO ATÉ A INTERSECÇÃO

                count += 1

            if self.ideal != 'IDEAL':
                if stage > 500:
                    self.show_popup('Erro. Ajustar valores!')
                    return

        # -----
        show_stage_num(dados)

        partial_stage = (xB - dados[-1][0]) / (dados[-3][0] - dados[-1][0])

        # Linhas de suporte para interpretação
        self.canvas.axes.plot([xD, xD], [0, xD], color='g', alpha=0.6, linestyle=':')
        self.canvas.axes.plot([xF, xF], [0, xF], color='g', alpha=0.6, linestyle=':')
        self.canvas.axes.plot([xB, xB], [0, xB], color='g', alpha=0.6, linestyle=':')

        # Config gráfico
        self.canvas.axes.minorticks_on()
        self.canvas.axes.grid(which='both', linestyle=':', alpha=0.6)
        self.canvas.axes.grid(visible=True, which='minor', linestyle=':', alpha=0.3)
        self.ui.horizontalLayout.addWidget(self.canvas)

        # Resultados
        res_stage = f'{str(stage)} ({str(round(stage - partial_stage, 2))})'
        self.ui.res_estagios.setText(res_stage)
        self.ui.res_alimentacao.setText(str(feed_stage))

        self.canvas.axes.margins(y=0.01)

        # self.canvas.axes.set_aspect('equal', adjustable='box')

        plt.close(fig)

    # /\/\ --------------------------   MCCABE FIM  --------------------------------- /\/\

    # ------------------ FUG INÍCIO ------------------

    def update_table(self):
        # Atualiza a listagem dos número no combobox, a partir do spinbox

        self.ui.comps_num_qBox.clear()

        # Após limpar, seleciona os valores do spin, e passa para uma lista -> combobox
        qtd = self.ui.qtd_comps.value()
        nums_list = [str(x) for x in range(1, qtd + 1)]
        self.ui.comps_num_qBox.addItems(nums_list)

        # Se a qtd comps selecionados > qtd números, retira os últimos componentes add
        if len(self.comps) > len(nums_list):
            aux = len(self.comps) - len(nums_list)
            for x in range(aux):
                self.comps.popitem()

        # Se a houver mais dados de frac molar do que a qtd de itens, retira os valores
        if len(self.mole_frac_vls) > qtd:
            for num in reversed(self.mole_frac_vls.keys()):
                del (self.mole_frac_vls[num])

                if len(self.mole_frac_vls) == qtd:
                    break

    def comp_selec(self):
        # Conecta os comps com seus índices
        index = self.ui.comps_num_qBox.currentIndex() + 1
        comp = self.ui.comps_qBox.currentText()

        self.comps[index] = comp

    def key_pos(self):
        # Cruza indíce com LK e HK; verifica duplicidade
        index = self.ui.comps_num_qBox.currentIndex() + 1

        if self.ui.hk_rdbtn.isChecked():
            self.key_comps['HK'] = index
            self.check_dupl('HK')
        elif self.ui.lk_rdbtn.isChecked():
            self.key_comps['LK'] = index
            self.check_dupl('LK')

    def check_dupl(self, key):
        # Confere e deleta se houver + de 1 LK/HK
        if len(self.key_comps) > 1:
            if self.key_comps['HK'] == self.key_comps['LK']:
                if key == 'HK':
                    del(self.key_comps['LK'])
                elif key == 'LK':
                    del(self.key_comps['HK'])

    def cnum_click(self):
        # Faz ajuste dos botões LK/HK, dependendo do num selecionado

        num = self.ui.comps_num_qBox.currentIndex() + 1

        # Se o número for LK/HK deixa marcado o radiobutton
        if num in list(self.key_comps.values()):
            if num == self.key_comps['LK']:
                self.ui.lk_rdbtn.setChecked(True)
            elif num == self.key_comps['HK']:
                self.ui.hk_rdbtn.setChecked(True)

        # Se não, desmarca todos
        else:
            self.ui.buttonGroup_3.setExclusive(False)
            self.ui.lk_rdbtn.setChecked(False)
            self.ui.hk_rdbtn.setChecked(False)
            self.ui.buttonGroup_3.setExclusive(True)

        # Mostra a fração dos componentes já indexados, ou 0 se ainda não tiver sido informado
        if num in self.mole_frac_vls.keys():
            self.ui.mole_frac.setValue(self.mole_frac_vls[num])
        else:
            self.ui.mole_frac.setValue(0)

        # Ajuste de nomes dos componentes selecionados
        if num in self.comps.keys():
            self.ui.comps_qBox.setCurrentText(self.comps[num])  # Nome
        else:
            self.ui.comps_qBox.setCurrentIndex(-1)             # Em branco

    def check_vls(self):
        # Faz a preview dos comps, chave e índice
        text = 'Componentes: '

        # Para cada índice (x[0])/comp no dicionário, add a string text
        if len(self.comps) <= self.ui.qtd_comps.value():
            for x in self.comps.items():
                try:
                    text += f'\n {x[0]}: {x[1]}: {self.mole_frac_vls[x[0]]}'
                except:
                    return

                # Se for comp chave, add LK/HK ao lado
                if x[0] in list(self.key_comps.values()):

                    # Pelo índice, pega a pos do comp chave e escreve HK/LK (key)
                    aux = list(self.key_comps.keys())[list(self.key_comps.values()).index(x[0])]
                    text += f' - {aux}'

            text += f'\nTotal frações: {round(sum(self.mole_frac_vls.values()), 2)}'

            self.ui.results_fug.setText(text)

    def f_molar(self):
        # Guarda a fração molar do comp no dict
        value = round(self.ui.mole_frac.value(), 2)
        comp = self.ui.comps_num_qBox.currentIndex() + 1

        if value != 0:
            self.mole_frac_vls[comp] = value

    def exec_fug(self):

        # Análise valores de entradas
        F = self.check_input(self.ui.F.text())
        ref_ratio = self.check_input(self.ui.ref_ratio.text())
        q = self.check_input(self.ui.q_fug.text())
        Tbolha = self.check_input(self.ui.Tbolha.text())
        Torvalho = self.check_input(self.ui.Torvalho.text())
        f_lk = round(self.ui.f_lk.value(), 2)
        f_hk = round(self.ui.f_hk.value(), 2)

        qtd = self.ui.qtd_comps.value()
        selec_comps = list(self.comps.values())

        try:
            pos_lk = self.key_comps['LK']
            pos_hk = self.key_comps['HK']
        except:
            self.show_popup('Falta determinar HK e/ou LK')
            return

        verificar = [F, ref_ratio, q, Tbolha, Torvalho, f_lk, f_hk]

        if (None in verificar) or (0 in verificar):
            self.show_popup('Faltam valores preenchidos')
            return
        elif len(selec_comps) != qtd:
            # Informados menos componentes do que o número informado
            self.show_popup('Número de componentes errado')
        elif (0 in self.mole_frac_vls) or (len(self.mole_frac_vls) < qtd):
            # Fração molar zerada ou faltando
            self.show_popup('Erro de valores')
        else:

            # Se a qtd de frações difere do n° de comps escolhidos ou a soma difere de 1, retorna
            if len(self.mole_frac_vls) != len(self.ui.comps_num_qBox) or round(sum(self.mole_frac_vls.values()), 2) != 1:
                return

            # Se LK e HK não são vizinhos, retorna
            if pos_lk != (pos_hk - 1):
                return

            self.calc_fug(F, ref_ratio, q, Tbolha, Torvalho, f_lk, f_hk, qtd, selec_comps, pos_lk, pos_hk)

    def calc_fug(self, F, ref_ratio, q, Tbolha, Torvalho, f_lk, f_hk, qtd, selec_comps, pos_lk, pos_hk):
        # Realiza os cálculos

        # Filtra o df pelos componentes escolhidos e formata ele
        df = self.dados_antoine.loc[selec_comps]
        df.reset_index(inplace=True)
        df['Index'] = np.arange(1, len(df) + 1)
        df.index = np.arange(1, len(df) + 1)

        # Add fração molar ao df
        df['z'] = df['Index'].map(self.mole_frac_vls)

        # Ajusta a qtd de itens do df e a atual informada
        if len(df) > qtd:
            diff = qtd - len(df)
            del_items = list(df.index)[-diff:]

            df.drop(del_items)
            df.reset_index()

        # Componentes chave e fração que sobra no topo e fundo
        df.at[pos_lk, 'xDD'] = F * df.at[pos_lk, 'z'] * f_lk  # LK que sobra no topo
        df.at[pos_lk, 'xBB'] = F * df.at[pos_lk, 'z'] * (1 - f_lk)

        df.at[pos_hk, 'xDD'] = F * df.at[pos_hk, 'z'] * (1 - f_hk)  # Chave pesada que sobra no topo da coluna
        df.at[pos_hk, 'xBB'] = F * df.at[pos_hk, 'z'] * f_hk

        for comp in df.index:
            if comp < pos_lk:
                # Comps + leves que LK zerados no fundo e determinação no topo
                df.at[comp, 'xDD'] = F * df.at[comp, 'z']
                df.at[comp, 'xBB'] = 0.0
            if comp > pos_hk:
                # Comps + pesados que HK zerados no topo e determinação no fundo
                df.at[comp, 'xDD'] = 0.0
                df.at[comp, 'xBB'] = F * df.at[comp, 'z']

        # Fração de destialdo e fundo de cada comp, após distr
        df['xD'] = df['xDD'] / sum(df['xDD'])
        df['xB'] = df['xBB'] / sum(df['xBB'])

        # Vazões fundo e topo
        B = sum(df['xBB'])
        D = sum(df['xDD'])

        # Antoine - Eq. (10)
        df['bolha'] = 10 ** (df['A'] - (df['B'] / (df['C'] + Tbolha)))      # Fundo
        df['orvalho'] = 10 ** (df['A'] - (df['B'] / (df['C'] + Torvalho)))  # Topo

        # Eq. (9) - Volatilidades relativas fundo e topo
        df['al_b'] = (df['bolha']) / (df['bolha'][pos_hk])
        df['al_d'] = (df['orvalho']) / (df['orvalho'][pos_hk])

        # Eq. (12) - volatilidade relativa média
        df['vol_rel_med'] = (df['al_d'] * df['al_b']) ** 0.5

        # Fenske - Eq. (11)
        dest = df['xD'][pos_lk] / df['xD'][pos_hk]
        fundo = df['xB'][pos_hk] / df['xB'][pos_lk]

        Nmin = (math.log(dest * fundo) / math.log(df['vol_rel_med'][pos_lk])) - 1

        # !! COMPONENTES DISTRIBUÍDOS
        new_df = df.copy()
        for row_num in new_df.index:
            if row_num not in [pos_lk, pos_hk]:

                # Eq. (13) - Componentes não chave - distribuição
                xidd_xibb = (df['vol_rel_med'][row_num] ** Nmin) * (df['xDD'][pos_hk] / df['xBB'][pos_hk])
                xibb = (F * df['z'][row_num]) / (xidd_xibb + 1)
                xidd = (F * df['z'][row_num] - xibb)

                new_df.at[row_num, 'xBB'] = round(xibb, 2)
                new_df.at[row_num, 'xDD'] = round(xidd, 2)

        # Atualização frações
        new_df['xD'] = new_df['xDD'] / sum(new_df['xDD'])
        new_df['xB'] = new_df['xBB'] / sum(new_df['xBB'])

        # UNDERWOOD
        calculo_theta = ''

        # Eq. (14) - Determinação theta
        for x in new_df.index:
            volatilidade = new_df['vol_rel_med'][x]
            fracao = new_df['z'][x]

            exp = f'({volatilidade}*{fracao})/({volatilidade}-y)'
            calculo_theta += exp

            if x < len(new_df.index):
                calculo_theta += '+'

        calculo_theta = f'{calculo_theta},1-{q}'
        y = sp.symbols('y')  # Theta

        calculo_theta_conv = sympify(f"Eq({calculo_theta})")
        results_theta = list(sp.solveset(calculo_theta_conv))

        # Validação valores theta
        for num in results_theta:
            if new_df['vol_rel_med'][pos_lk] > num > 1:
                theta = float(num)  # Valor correto

        # Eq. (15) - Rmin: Underwood
        calculo_rmin = ''

        for x in new_df.index:
            volatilidade = new_df['vol_rel_med'][x]
            fracao = new_df['xD'][x]

            exp2 = f'({volatilidade}*{fracao})/({volatilidade}-{theta})'
            calculo_rmin += exp2

            if x < len(new_df.index):
                calculo_rmin += '+'

        calculo_rmin = f'{calculo_rmin},1+R'
        R = sp.symbols('R')

        calculo_rmin_conv = sympify(f"Eq({calculo_rmin})")
        results_rmin = sp.solveset(calculo_rmin_conv)
        Rmin = float(list(results_rmin)[0])  # Valor certo

        # Gilliland
        R = Rmin * ref_ratio
        N = sp.symbols('N')

        # Eq. 16 e 17
        x_eq17 = (R - Rmin) / (R + 1)
        gill = sp.Eq(1-math.exp(((1 + 54.4 * x_eq17)/(11 + 117.2 * x_eq17)) * ((x_eq17 - 1) / x_eq17 ** 0.5)),
                     (N - Nmin) / (N + 1))
        solution_gilliland = sp.solveset(gill)
        N = float(list(solution_gilliland)[0])

        # Kirkbride - Eq. (18)
        nrns = ((new_df['z'][pos_hk] / new_df['z'][pos_lk]) * ((new_df['xB'][pos_lk] / new_df['xD'][pos_hk]) ** 2) * (B / D)) ** 0.206
        Nr = (nrns * N) / (1 + nrns)

        # RESULTADOS
        final_results = f'Nmin: {round(Nmin, 3)}\nRmin: {round(Rmin, 3)}\nN: {round(N, 3)}\nNr: {round(Nr, 3)}\n'
        final_results += '-' * 30
        final_results += '\n' + str(new_df[['Substance', 'xDD', 'xBB', 'xD', 'xB']].round(4)) + '\n'
        final_results += '-' * 30
        final_results += '\nLEGENDA:'
        final_results += '\nNmin: número de estágios mínimos. Rmin: razão de refluxo mínima. N: número de estágios. ' \
                         'Nr: número de estágios de retificação (local de alimentação)'

        self.ui.results_fug.clear()
        self.ui.results_fug.setText(final_results)

    # ------------------- FUG FIM --------------------


if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_win = MainWindow()
    main_win.show()
    sys.exit(app.exec_())
