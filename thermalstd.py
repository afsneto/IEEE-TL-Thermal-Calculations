# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:58:03 2018

@author: afsn3
"""

import math as mt
import numpy as np
import pandas as pd
import cmath as cm
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from cable import loadcable
from tqdm import tqdm
# from cable import optimalcable


class climavars:

    def __init__(self, vw, em, ab, tamb, zl, lat, atm, he, phi, hour, nday):

        self.vw = vw  # Wind speed (m/s)
        self.em = em  # Emissivity
        self.ab = ab  # Solar absorptivity
        self.tamb = tamb  # Ambient air temperature (°C)
        self.zl = zl  # Azimuth of line
        self.lat = lat  # Latitude - positive values North, negative South
        self.atm = atm  # 1 to Clear atmosphere and 2 to Industrial atmosphere
        self.he = he  # Average conductor elevation (m)
        self.phi = phi  # Angle between the wind direction and the conductor  axis
        self.hour = hour  # hour day
        self.nday = nday  # day of the year


class cablevars:
    def __init__(self, db_cable, cable_type, cable_name):
        self.db_cable = db_cable
        self.cable_type = cable_type
        self.cable_name = cable_name


class Std7382006:
    def __init__(self, climavars, cablevars):

        self.vw = climavars.vw  # Wind speed (m/s)
        self.em = climavars.em  # Emissivity
        self.ab = climavars.ab  # Solar absorptivity
        self.tamb = climavars.tamb  # Ambient air temperature (°C)
        self.zl = climavars.zl  # Azimuth of line
        self.lat = climavars.lat  # Latitude - positive values North, negative South
        self.atm = climavars.atm  # 1 to Clear atmosphere and 2 to Industrial atmosphere
        self.he = climavars.he  # Average conductor elevation (m)
        self.phi = climavars.phi  # Angle between the wind direction and the conductor  axis
        self.hour = climavars.hour  # hour day
        self.nday = climavars.nday  # day of the year

        self.db_cable = cablevars.db_cable
        self.cable_type = cablevars.cable_type
        self.cable_name = cablevars.cable_name

        self.paramcable = self.selectcable()

    def selectcable(self):
        # Loads the cable data to calculate the long and short-term current.
        cable_data = loadcable(self.db_cable).filtercable(
            self.cable_type, self.cable_name)
        return cable_data

    def graphcable(self):
        x, y = self._pointscable()

        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(x, y)
        ax.set_xlabel('Current (A)')
        ax.set_ylabel('Temperature (°C)')
        plt.grid()
        plt.show()

    def _resistances_cable(self):

        if self.paramcable['ElectRes_CC_20'].values[0].any():
            # km to meters /1000
            r1 = self.paramcable['ElectRes_CC_20'].values[0] / 1000
            t1 = 25
        else:
            r1 = self.paramcable['ElectRes_CA_25'].values[0] / 1000
            t1 = 20

        if self.paramcable['ElectRes_CA_75'].values[0].any():
            r2 = self.paramcable['ElectRes_CA_75'].values[0] / 1000
            t2 = 50
        else:
            r2 = self.paramcable['ElectRes_CA_50'].values[0] / 1000
            t2 = 75

        return r1, t1, r2, t2

    def _qc(self, tc):
        """
        Calculating Convection heat loss (qc)
        :return: qc
        """
        if self.atm == 1:
            self.atmType = 'Clear'
        else:
            self.atmType = 'Industrial'
        '''
        Natural convection heat loss (qcn)
        With zero wind speed
        pf =  Air density (kg/m³)
        '''

        # For the cases tc - tamb < 0, qcn = 0
        # There is not natural convection

        t_film = (tc + self.tamb)/2

        pf = (1.293 - 1.525e-4 * self.he + 6.379e-9 *
              self.he ** 2) / (1 + 0.00367 * t_film)

        # For the cases tc - tamb < 0, qcn = 0
        # There is not natural convection
        # del_t = (tc - self.tamb).clip(0)

        delt_t = tc - self.tamb
        if delt_t < 0:
            delt_t = 0
        else:
            delt_t = delt_t

        delt_t = tc - self.tamb

        def powernumpy(x, y):
            return np.sign(x) * (np.abs(x)) ** y

        diamcable = self.paramcable['Diam_Total'].values[0]
        # Deal with RuntimeWarning: invalid value encountered in power
        # qcn = 0.02050 * pf ** 0.5 * \
        #     self.paramcable['Diam_Total'].values[0] ** 0.75 * delt_t ** 1.25
        qcn = 0.02050 * \
            powernumpy(pf, 0.5) * powernumpy(diamcable, 0.75) * \
            powernumpy(delt_t, 1.25)

        '''
        Forced convection heat loss (qc1 e qc2)
        qc1 = Applies at low winds
        qc2 = Applies at high winds
        The larger of the two calculated convection heat loss rates is used
        mf = Dynamic viscosity of air (Pa-s)
        kf = Thermal conductivity of air
        '''
        mf = (1.458e-6 * (t_film + 273) ** 1.5) / (t_film + 383.4)

        kf = 2.424e-2 + 7.477e-5 * t_film - 4.407e-9 * t_film ** 2

        '''
        k_angle = Wind direction factor
        phi = Angle between the wind direction and the conductor axis
        '''
        phirad = mt.radians(self.phi)
        k_angle = 1.194 - mt.cos(phirad + 0.194 *
                                 mt.cos(2 * phirad) + 0.368 * mt.sin(2 * phirad))

        qc1 = 1.01 + 0.0372 * ((diamcable * pf * self.vw) /
                               mf) ** 0.52 * kf * k_angle * delt_t

        qc2 = 0.0119 * ((diamcable * pf * self.vw) /
                        mf) ** 0.6 * kf * k_angle * delt_t

        if self.vw == 0:
            qc = qcn
        else:
            qc = np.maximum(qc1, qc2)

        return qc

    def _qr(self, tc):
        """
        Calculate radiated heat loss (qr)
        :param self:
        :return: qr
        """
        qr = 0.0178 * self.paramcable['Diam_Total'].values[0] * self.em * (
            (((tc + 273) / 100) ** 4) - ((self.tamb + 273) / 100) ** 4)

        return qr

    def _qs(self, n_day=161):
        """
        Calculate the solar heat gain of a conductor.
        :param n_day: Number of the year
        :return: qs: Solar heat gain
        """

        if n_day == '':
            n_day = self.nday
        else:
            n_day = n_day

        w_hour = (self.hour - 12) * 15

        # Include an option to choose between manual day put or calculating for an
        # annual peak solar heat input

        # Use a function that find a value of N to evaluating the large value of gama
        # sin (90°) = 1

        # test for the June 10
        # gama = 23.458 * sin(((284 + N)/365)*360)
        # June 10 = 161° day year

        gama = 23.458 * mt.sin((mt.radians(284 + n_day) / 365) * 360)

        # Hc = Altitude of the sun
        hc_rad = mt.asin(mt.cos(mt.radians(self.lat)) * mt.cos(mt.radians(gama)) * mt.cos(
            mt.radians(w_hour)) + mt.sin(mt.radians(self.lat)) * mt.sin(mt.radians(gama)))
        hc = mt.degrees(hc_rad)

        # sav = Solar azimuth variable
        sav = mt.sin(mt.radians(w_hour)) / (mt.sin(mt.radians(self.lat)) * mt.cos(mt.radians(w_hour))
                                            - mt.cos(mt.radians(self.lat)) * mt.tan(mt.radians(gama)))

        # sac = Solar azimuth constant
        sac = None
        if -180 <= w_hour < 0:
            if sav >= 0:
                sac = 0
            else:
                sac = 180
        elif 0 <= w_hour <= 180:
            if sav >= 0:
                sac = 180
            else:
                sac = 360

        # Zc = solar azimuth angle
        zc = sac + 180 / mt.pi * (mt.atan(sav))

        # Solar heat gain
        theta = mt.acos(mt.cos(hc_rad) *
                        mt.cos(mt.radians(zc) - mt.radians(self.zl)))

        # Total heat flux elevation correction factor

        qs = None

        if self.atm == 1:
            a = -42.2391
            b = 63.8044
            c = -1.9220
            d = 3.46921e-2
            e = -3.61118e-4
            f = 1.94318e-6
            g = -4.07608e-9

            qs_total_solar = a + b * hc + c * hc ** 2 + d * \
                hc ** 3 + e * hc ** 4 + f * hc ** 5 + g * hc ** 6
            k_solar = 1 + 1.148e-4 * self.he - 1.108e-8 * self.he ** 2
            qse = k_solar * qs_total_solar
            qs = self.ab * qse * mt.sin(theta) * \
                self.paramcable['Diam_Total'].values[0] / 1000

        elif self.atm == 2:
            a = 53.1821
            b = 14.2110
            c = 6.6138e-1
            d = -3.1658e-2
            e = 5.4654e-4
            f = -4.3446e-6
            g = 1.3236e-8

            qs_total_solar = a + b * hc + c * hc ** 2 + d * \
                hc ** 3 + e * hc ** 4 + f * hc ** 5 + g * hc ** 6
            k_solar = 1 + 1.148e-4 * self.he - 1.108e-8 * self.he ** 2
            qse = k_solar * qs_total_solar
            qs = self.ab * qse * mt.sin(theta) * \
                self.paramcable['Diam_Total'].values[0] / 1000
        return qs

    def _pointscable(self):
        tempinterval = np.arange(10, 150, 1)
        x = [self.current(i) for i in tempinterval]
        return x, tempinterval

    def current(self, tc):

        qc = self._qc(tc)
        qr = self._qr(tc)
        qs = self._qs()

        r1, t1, r2, t2 = self._resistances_cable()
        rf = r1 + ((r2 - r1) / (t2 - t1)) * (tc - t1)
        final_current = cm.sqrt((qc + qr - qs) / rf)
        return np.real(final_current)

    def temp(self, current, pointslist):
        if pointslist:
            x, y = pointslist
        else:
            x, y = self._pointscable()
        # POLINOMIO 4 GRAU MELHOR COMPORTAMENTO COM CURVA
        pol = np.polyfit(x, y, 4)
        temp = np.polyval(pol, current)
        return temp

    def new_resistance(self, tf=50):
        if self.paramcable['ElectRes_CC_20'].values[0].any():
            r1 = self.paramcable['ElectRes_CC_20'].values[0] / 1000
            t1 = 20
            rf = r1 * (1 + (self._calculate_alpha() * (tf - t1)))
        else:
            r1, t1, r2, t2 = self._resistances_cable()
            rf = r1 + ((r2 - r1) / (t2 - t1)) * (tf - t1)
        return rf  # ohm / m

    def power_loss(self, rf, current=0):
        # final_temp = self.calculate_temp(current)
        # rtc_1, tc_1, rtc_2, tc_2 = self._resistances_cable()
        if rf:
            rf = rf
        else:
            rf = self.new_resistance()
        # rf em metros * 1e3 = rf em km
        # power_loss em W / km , divido por 1e6 = MW / km
        power_loss = 3 * (current ** 2 * (rf * 1e3) / (1e6))

        return power_loss

    def _cable_description(self):
        if self.paramcable.Name.isnull().any():
            cable_description = (str(self.paramcable['Type'].values[0]) + ' ' + str(self.paramcable['AWG'].values[0]) + ' MCM ' +
                                 str(self.paramcable['Diam_Total'].values[0]) + ' mm')
        elif self.paramcable.AWG.isnull().any():
            cable_description = (str(self.paramcable['Type'].values[0]) + ' ' + str(self.paramcable['Name'].values[0]) + ' MCM ' +
                                 str(self.paramcable['Diam_Total'].values[0]) + ' mm')
        else:
            if self.paramcable['Name'].values[0] == self.paramcable['AWG'].values[0]:
                cable_description = (str(self.paramcable['Type'].values[0]) + ' ' + str(self.paramcable['AWG'].values[0]) + ' MCM ' +
                                     str(self.paramcable['Diam_Total'].values[0]) + ' mm')
            else:
                cable_description = (str(self.paramcable['Type'].values[0]) + ' ' + str(self.paramcable['Name'].values[0]) + ' ' +
                                     str(self.paramcable['AWG'].values[0]) + ' MCM ' +
                                     str(self.paramcable['Diam_Total'].values[0]) + ' mm')
        return cable_description

    def _calculate_alpha(self):
        cables_type = self.paramcable['Type'].values[0]
        if cables_type == 'ACSR':
            alpha = 0.00403
        elif cables_type == 'AAAC_6201':
            alpha = 0.00347
        elif cables_type == 'AAAC_1120':
            alpha = 0.0039
        else:
            alpha = 0.00403  # ACAR and AAC
        return alpha


class losspower(Std7382006):

    def __init__(self, dfproduction, voltage, powerfactor, climavars, cablevars, outnetlimit, extline):
        self.df = dfproduction
        self.voltage = voltage
        self.powerfactor = powerfactor
        self.outnetlimit = outnetlimit
        self.extline = extline
        super().__init__(climavars, cablevars)

    def _current(self, power):
        return (power * 1e6 * self.powerfactor) / (self.voltage * mt.sqrt(3))

    def _completedf(self, df):
        # df['I(A)'] = np.zeros(df.shape[0])
        # df['Cable Temp'] = np.zeros(df.shape[0])
        # df['Joule Loss'] = np.zeros(df.shape[0])
        # df['Out Net Prod'] = np.zeros(df.shape[0])
        # df['% JL'] = np.zeros(df.shape[0])

        newr = [self.new_resistance(i) for i in df.loc[:, 'Temperature']]

        df['I(A)'] = [self._current(i) for i in df.loc[:, 'Net Prod']]
        pointslist = self._pointscable()
        df['Cable Temp'] = [self.temp(i, pointslist)
                            for i in df.loc[:, 'I(A)']]
        df['Joule Loss'] = [self.power_loss(i, j) * self.extline
                            for i, j in zip(newr, df.loc[:, 'I(A)'])]
        df['Out Net Prod'] = df['Net Prod'] - df['Joule Loss']
        df['% JL'] = np.nan_to_num(df['Joule Loss']/df['Out Net Prod']) * 100
        return df

    def initdf(self):
        df = self.df
        newdf = self._completedf(df)
        return newdf

    def opt_netprod(self, netprod):

        obj = self._calculationvars(netprod)
        if obj <= 0:
            pass
        else:
            while obj > 0:
                obj = self._calculationvars(netprod)
                netprod -= 0.05

        return netprod

    def finaldf(self):

        df = self.initdf()
        df['Net Prod'] = np.clip(df['Net Prod'], 0, 62)
        # df = self.initdf()
        # vfunc = np.vectorize(self.opt_netprod)
        # df['Net Prod'] = vfunc(df['Net Prod'])
        # df['Net Prod'] = [self.opt_netprod(x)
        #                   for index, x in tqdm(np.ndenumerate(df['Net Prod']))]
        return self._completedf(df).to_excel('dfteste.xlsx')

    def _calculationvars(self, netprod):
        if netprod > self.outnetlimit:

            pointslist = self._pointscable()
            new_cur = self._current(netprod)
            new_cabletemp = self.temp(new_cur, pointslist)
            new_r = self.new_resistance(new_cabletemp) * 1000  # loss / km
            new_jouleloss = self.power_loss(new_r, new_cur) * self.extline
            new_outnet = netprod - new_jouleloss

            obj = (self.outnetlimit - new_outnet) ** 2

        else:
            obj = 0

        return obj

    def initdf_tqdm(self):
        df = self.df

        df['I(A)'] = np.zeros(df.shape[0])
        df['Cable Temp'] = np.zeros(df.shape[0])
        df['Joule Loss'] = np.zeros(df.shape[0])
        df['Out Net Prod'] = np.zeros(df.shape[0])
        df['% JL'] = np.zeros(df.shape[0])

        newr = []
        pointslist = self._pointscable()

        for i in tqdm(range(df.shape[0])):
            newr.append(self.new_resistance(
                df.loc[df.index[i], ['Temperature']].values[0]))

            df.loc[df.index[i], ['I(A)']] = self._current(
                df.loc[df.index[i], ['Net Prod']])

            df.loc[df.index[i], ['Cable Temp']] = self.temp(
                df.loc[df.index[i], ['I(A)']], pointslist)

            df.loc[df.index[i], ['Joule Loss']] = self.power_loss(
                newr[i], df.loc[df.index[i], ['I(A)']])

            df.loc[df.index[i], ['Out Net Prod']] = df.loc[df.index[i],
                                                           ['Net Prod']] + df.loc[df.index[i], ['Joule Loss']]

            df.loc[df.index[i], ['% JL']] = df.loc[df.index[i], [
                'Joule Loss']] / df.loc[df.index[i], ['Out Net Prod']]

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
