import math as mt
import numpy as np
import pandas as pd
from tqdm import tqdm
import swifter

from thermalstd import Std7382006

# TODO COMPLETE DATAFRAME GROSS AND NET POWER WITH CURRENT, CABLE TEMP, JOULE LOSS AND LIMITED POWER


class losspower(Std7382006):

    def __init__(self, dfproduction, voltage, powerfactor):
        self.df = dfproduction
        self.voltage = voltage
        self.powerfactor = powerfactor

    def _current(self, power):
        return (power * 1e6 * self.powerfactor) / (self.voltage * mt.sqrt(3))

    def completedf(self):
        df = self.df

        df['I(A)'] = np.zeros(df.shape[0])
        df['Cable Temp'] = np.zeros(df.shape[0])
        df['Joule Loss'] = np.zeros(df.shape[0])
        df['Out Net Prod'] = np.zeros(df.shape[0])

        tqdm.pandas()
        df.loc[:, 'I(A)'] = df['Net Prod'].swifter.apply(self._current)
        df.loc[:, 'Cable Temp'] = df['I(A)'].swifter.apply(self.temp)
        return df

    # def cablechoice(self, cable_type, cable_name, voltage, extline):
    #         # TODO
    #         # if csvfile does not exist, loadpandasfile and execute cablechoicew
    #         # else
    #         # df = pd.read_csv(csvfile)

    #         df = self.df
    #         self.cable_type = cable_type
    #         self.cable_name = cable_name
    #         self.voltage = voltage
    #         self.extline = extline
    #         filename_csv = cable_type + '_' + str(cable_name) + '.csv'
    #         filename_excel = cable_type + '_' + str(cable_name) + '.xlsx'

    #         data = Std7382006(cable_type, cable_name)

    #         # Calculating the current for each net production
    #         # CALCULAR CORRENTE PARA df['Out Power Gross']*1e6/(voltage * np.sqrt(3))
    #         df['I(A)'] = df['Net Prod']*1e6/(voltage * np.sqrt(3))
    #         df['Cable Temp'] = np.zeros(df.shape[0])

    #         df['Joule Loss'] = np.zeros(df.shape[0])

    #         df['Out Net Prod'] = np.zeros(df.shape[0])
    #         df['Joule Loss %'] = np.zeros(df.shape[0])

    #         my_list = range(df.shape[0])

    #         print('Processing {}'.format(filename_csv))

    #         for i in tqdm(my_list):
    #             # df.iloc[i, 4]  # Wind speed
    #             # df.iloc[i, 2]  # Temperature
    #             data.climadata(df.iloc[i, 4],
    #                            self.em,
    #                            self.ab,
    #                            df.iloc[i, 2],
    #                            self.zl,
    #                            self.lat,
    #                            self.atm,
    #                            self.he,
    #                            self.phi,
    #                            self.hour,
    #                            self.nday)
    #             # Cable Temp
    #             df.iloc[i, 8] = data.calculate_temp(df.iloc[i, 7])

    #             # Joule Loss (MW)
    #             cable_res = data.calculate_new_resistance(df.iloc[i, 8], mode='ac')
    #             if df.iloc[i, 6] == 0:
    #                 df.iloc[i, 9] = 0
    #             else:
    #                 df.iloc[i, 9] = 3 * ((cable_res * 1000) * df.iloc[i, 7]
    #                                      ** 2) * extline / 1e6  # MW

    #             # Out Net Prod = Net Prod + Joule Loss (MW)
    #             if df.iloc[i, 6] > 60:
    #                 df.iloc[i, 10] = 60 + df.iloc[i, 9]
    #             else:
    #                 df.iloc[i, 10] = df.iloc[i, 6] + df.iloc[i, 9]

    #             # Joule Loss / Out Net Prod
    #             if df.iloc[i, 9] == 0:
    #                 df.iloc[i, 11] = 0
    #             else:
    #                 df.iloc[i, 11] = df.iloc[i, 9] / (df.iloc[i, 10]) * 100

    #         df.to_csv(filename_csv)
    #         print('File {} saved'.format(filename_csv))
    #         # TODO RESOLVE TypeError trying save .xlsx file
    #         # TypeError: got invalid input value of type <class 'xml.etree.ElementTree.Element'>, expected string or Element
    #         # update openpyxl from 3.0.2 to 3.0.3
    #         # pip install openpyxl==3.0.3
    #         df.to_excel(filename_excel)
    #         print('File {} saved\n'.format(filename_excel))

    #         # print ouput
    #         rtc_1, tc_1, rtc_2, tc_2 = data._resistances_cable()

    #         print('CABLE CHARACTERISTICS')
    #         print('TYPE: {}\nNAME: {}'.format(cable_type, cable_name))
    #         print('RESISTANCE {0:.1f}: {1:.4f} ohm/km'.format(tc_1, rtc_1 * 1000))
    #         print('RESISTANCE {0:.1f}: {1:.4f} ohm/km'.format(tc_2, rtc_2 * 1000))
    #         print('\n')

    #     def cableverify(self, csvfile, longtermtemp, lossper):
    #         '''
    #         VERIFICAMOS SE A MÉDIA DAS PERDAS ELÉTRICAS CA ESTÁ ABAIXO DO VALOR SOLICITADO NO CONTRATO

    #         CORRENTE MÁXIMA ENCONTRADA 666 A = CORRENTE DE CURTA DURAÇÃO
    #         PARA UMA TEMPERATURA DE LONGA DURAÇÃO IGUAL A 60°C, CONFORME REN 191
    #         TEMOS UM FATOR K IGUAL A 1.26, O QUE RESULTA EM UMA CORRENTE DE LONGA DURAÇÃO = 529 A
    #         VERIFICAMOS PARA QUAL PORCENTAGEM (QUANTILE) RESULTA UMA CORRENTE DE 529 A AO lONGO DOS X ANOS DE AMOSTRA
    #         VERIFICAMOS QUAL A TEMPERATURA ALCANÇADA PARA O CABO PARA AS CONDIÇÕES DE LONGA E CURTA DURAÇÃO SEGUNDO A IEEE
    #         COMPARAMOS COM OS VALORES ENCONTRADOS DE TEMPERATURA CONSIDERANDO AS VARIÁVEIS CLIMÁTICAS DA AMOSTRA
    #         '''
    #         df = pd.read_csv(csvfile)
    #         meanlossjoule = df['Joule Loss %'].mean()

    #         print('DESCRIPTION OF I(A) AND CABLE TEMP')
    #         print(df[['I(A)', 'Cable Temp', 'Joule Loss %']].describe())

    #         # REN 191
    #         x = [50, 55, 60, 64, 65, 70, 75, 80, 90]
    #         y = [1.42, 1.33, 1.26, 1.24, 1.23, 1.19, 1.17, 1.15, 1.12]
    #         # POLINOMIO 4 GRAU MELHOR COMPORTAMENTO COM CURVA
    #         pol = np.polyfit(x, y, 4)
    #         k_factor = np.polyval(pol, longtermtemp)

    #         print('\nCONDITION 1')
    #         if meanlossjoule <= lossper:
    #             print('OK. THE AC LOSS {0:.2f} % IS <= {1:.2f} %.'.format(
    #                 meanlossjoule, lossper))
    #         else:
    #             print('NG. THE AC LOSS {0:.2f} % IS >= {1:.2f} %.'.format(
    #                 meanlossjoule, lossper))

    #         # CALCULATING LONG AND SHORT-TERM CURRENTS DATABASE
    #         shortcurrent = df['I(A)'].max()
    #         longcurrent = shortcurrent / k_factor

    #         shortcurrentadopt = round((shortcurrent + 10) / 10.0) * 10
    #         longcurrentadopt = round((longcurrent + 10) / 10.0) * 10

    #         templongcurrent = df['Cable Temp'][df['I(A)'] >= longcurrent].values[0]
    #         tempshortcurrent = df['Cable Temp'][df['I(A)']
    #                                             == shortcurrent].values[0]

    #         # CALCULATING LONG AND SHORT-TERM CURRENTS STANDARD 738- 2006
    #         data = Std7382006(self.cable_type, self.cable_name)
    #         vw = 1
    #         tamb = round(df['Temperature'].mean())
    #         data.climadata(vw, self.em, self.ab, tamb,
    #                        self.zl, self.lat, self.atm, self.he, self.phi,
    #                        self.hour, self.nday)
    #         templongcurrentstd = data.calculate_temp(longcurrent)
    #         tempshortcurrentstd = data.calculate_temp(shortcurrent)

    #         # simpleplot
    #         data.simpleplot(longcurrentadopt)
    #         data.simpleplot(shortcurrentadopt)

    #         print('\nCONDITION 2')
    #         print('LONG-TERM CURRENT IN DATABASE: {0:.2f} A'.format(longcurrent))
    #         print('CABLE TEMPERATURE IN DATABASE: {0:.2f}°C'.format(
    #             templongcurrent))

    #         print('LONG-TERM CURRENT ADOPT: {0:.2f} A'.format(longcurrentadopt))
    #         print(
    #             'CABLE TEMPERATURE IN IEEE STD 738-2006: {0:.2f} °C'.format(templongcurrentstd))

    #         print(
    #             '\nSHORT-TERM CURRENT IN DATABASE: {0:.2f} A'.format(shortcurrent))
    #         print('CABLE TEMPERATURE: {0:.2f}°C'.format(tempshortcurrent))

    #         print('SHORT-TERM CURRENT ADOPT: {0:.2f} A'.format(shortcurrentadopt))
    #         print(
    #             'CABLE TEMPERATURE IN IEEE STD 738-2006: {0:.2f} °C'.format(tempshortcurrentstd))
    #         print('\n')

    #         print('DAY OF SHORT-TERM CURRENT')
    #         print(df.iloc[df['I(A)'].idxmax()])

    #         per = np.arange(0, 1.1, 0.1)
    #         quantiles = df[['I(A)', 'Cable Temp']].quantile(per)

    #         # TODO CALCULATE FOR WHICH CURRENT VALUE >= LONGTERMTEMP

    #         print('\nPERCENTILE STATS')
    #         print(quantiles)
