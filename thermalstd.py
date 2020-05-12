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
import seaborn as sns
from pathlib import Path
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


class linevars:
    def __init__(self, dfproduction, voltage, powerfactor, climavars, cablevars, extline, maxnetprod=None, outnetlimit=None):
        self.df = dfproduction
        self.voltage = voltage
        self.powerfactor = powerfactor
        self.outnetlimit = outnetlimit
        self.extline = extline
        self.maxnetprod = maxnetprod


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
        ax.set_title('{}_{}'.format(self.cable_type,
                                    self.cable_name))
        plt.grid()
        plt.savefig('{}_{}.png'.format(self.cable_type,
                                       self.cable_name))

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
        tempcable_interval = np.arange(10, 150, 1)
        x = [self.current(i) for i in tempcable_interval]
        return x, tempcable_interval

    def current(self, tc):

        qc = self._qc(tc)
        qr = self._qr(tc)
        qs = self._qs()

        r1, t1, r2, t2 = self._resistances_cable()
        rf = r1 + ((r2 - r1) / (t2 - t1)) * (tc - t1)
        final_current = cm.sqrt((qc + qr - qs) / rf)
        return np.real(final_current)

    def temp(self, current, pointslist=''):
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

    def __init__(self, climavars, cablevars, linevars):

        self.df = linevars.df
        self.voltage = linevars.voltage
        self.powerfactor = linevars.powerfactor
        self.outnetlimit = linevars.outnetlimit
        self.extline = linevars.extline
        self.maxnetprod = linevars.maxnetprod
        super().__init__(climavars, cablevars)

    def _current(self, power):
        return (power * 1e6) / (self.voltage * 1e3 * mt.sqrt(3) * self.powerfactor)

    def _print_info_cable(self):
        cable_data = self.selectcable()
        print('TYPE: {}'.format(cable_data['Type'].values[0]))
        print('NAME: {}'.format(cable_data['Name'].values[0]))
        print('AWG: {}'.format(cable_data['AWG'].values[0]))

    def _completedf(self, df):

        newr = [self.new_resistance(i) for i in df.loc[:, 'Temperature']]

        df['I(A)'] = [self._current(i) for i in df.loc[:, 'Net Prod']]
        pointslist = self._pointscable()
        df['Cable Temp'] = [self.temp(i, pointslist)
                            for i in df.loc[:, 'I(A)']]
        df['Joule Loss'] = [self.power_loss(i, j) * self.extline
                            for i, j in zip(newr, df.loc[:, 'I(A)'])]
        df['Out Net Prod'] = df['Net Prod'] - df['Joule Loss']
        df['% JL'] = np.nan_to_num(df['Joule Loss']/df['Net Prod']) * 100
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
                netprod -= 1

        return netprod

    def limit_maxnetprod(self):

        df = self.initdf()
        if self.maxnetprod:
            df['Net Prod'] = np.clip(df['Net Prod'], 0, self.maxnetprod)
        else:
            pass
        df = self._completedf(df)
        # df = self.initdf()
        # vfunc = np.vectorize(self.opt_netprod)
        # df['Net Prod'] = vfunc(df['Net Prod'])
        # df['Net Prod'] = [self.opt_netprod(x)
        #                   for index, x in tqdm(np.ndenumerate(df['Net Prod']))]
        # if outputxlsx == '':
        #     pass
        # else:
        #     df.to_excel(outputxlsx)

        return df

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


class analysis(losspower):

    def __init__(self, climavars, cablevars, linevars, savexlsx=False):
        # CABLE DESCRIPTION
        super().__init__(climavars, cablevars, linevars)
        rtc_1, tc_1, rtc_2, tc_2 = self._resistances_cable()

        self.df = self.limit_maxnetprod()
        self.cable_type = cablevars.cable_type
        self.cable_name = cablevars.cable_name
        self.linevars = linevars
        if savexlsx:
            # REMOVE TIMEZONE BEFORE EXPORT TO EXCEL
            filename = '{}_{}.xlsx'.format(
                cablevars.cable_type, cablevars.cable_name)
            self.df.index = self.df.index.tz_localize(None)
            self.df.to_excel(filename)

        print('\n')
        print('-' * 55)
        print('CABLE CHARACTERISTICS')
        self._print_info_cable()
        print('RESISTANCE {0:.1f}: {1:.4f} ohm/km'.format(tc_1, rtc_1 * 1000))
        print('RESISTANCE {0:.1f}: {1:.4f} ohm/km'.format(tc_2, rtc_2 * 1000))

    def conditions(self, lossper=1.3):
        if self.linevars.maxnetprod:
            self.cutconditions(lossper)
        else:
            self.nocutconditions(lossper)

    def cutconditions(self, lossper):  # TODO INSERT PARAMETER maxtemp (90 or 94 °C)
        # CONDIÇÃO 1 PERDAS MÉDIAS % JOULE ABAIXO DO VALOR MÁXIMO
        meanlossjoule = self.df['% JL'].mean()
        print('\n***CONDITION 1***\n')
        if meanlossjoule <= lossper:
            print('OK. THE AC LOSS MEAN {0:.2f} % IS <= {1:.2f} %.'.format(
                meanlossjoule, lossper))
        else:
            print('NG. THE AC LOSS MEAN {0:.2f} % IS >= {1:.2f} %.'.format(
                meanlossjoule, lossper))

        cur_longdur = self.df['I(A)'].max()

        # CALCULATING THE TEMPERATURE FOR LONG DURATION CURRENT

        # Fixing the ambient temperature to a constant value
        climavars.tamb = round(self.df['Temperature'].groupby(
            self.df.index.day).max().mean())

        print("Média máxima dia temperatura: {} °C".format(climavars.tamb))

        temp_longdur = self.temp(cur_longdur)

        # REN 191
        x = [50, 55, 60, 64, 65, 70, 75, 80, 90]
        y = [1.42, 1.33, 1.26, 1.24, 1.23, 1.19, 1.17, 1.15, 1.12]
        # POLINOMIO 4 GRAU MELHOR COMPORTAMENTO COM CURVA
        pol = np.polyfit(x, y, 4)
        k_factor = np.polyval(pol, temp_longdur)
        cur_shortdur = cur_longdur * k_factor
        temp_shortdur = self.temp(cur_shortdur)

        # CALCULATING LONG AND SHORT-TERM CURRENTS STANDARD 738- 2006
        print('\n***CONDITION 2***\n')

        print('LONG-TERM CURRENT: {0:.2f} A'.format(cur_longdur))
        print(
            'CABLE TEMPERATURE IN IEEE STD 738-2006: {0:.2f} °C'.format(temp_longdur))

        if (temp_longdur > 94):
            print('WARNING: CABLE TEMPERATURE ABOVE 94°C. NG')
        else:
            print('OK')

        # CALCULATING LONG AND SHORT-TERM CURRENTS DATABASE
        print('SHORT-TERM CURRENT: {0:.2f} A'.format(cur_shortdur))
        print(
            'CABLE TEMPERATURE IN IEEE STD 738-2006: {0:.2f} °C'.format(temp_shortdur))
        if (temp_shortdur > 94):
            print('WARNING: CABLE TEMPERATURE ABOVE 94°C. NG')
        else:
            print('OK')

        # CALCULO CORRENTE DE LONGA CONDIÇÃO 3
        print('\n***CONDITION 3***\n')

        print('MAXIMUM TRANSMISSION LINE LOSS POWER: {0:.2f} MW'.format(
            self.df['Joule Loss'].max()))
        print('PERCENTAGE OF MAXIMUM TRANSMISSION LINE LOSS POWER: {0:.2f} %'.format(
            self.df['% JL'].max()))
        print('PERCENTAGE OF MEAN TRANSMISSION LINE LOSS POWER: {0:.2f} %'.format(
            self.df['% JL'].mean()))
        print('MAXIMUM NET POWER FROM POWERPLANT: {0:.2f} MW'.format(
            self.df['Net Prod'].max()))
        print('MAXIMUM POWER DELIVERED TO SUBSTATION: {0:.2f} MW'.format(
            self.df['Out Net Prod'].max()))

        print('\f', end='')

    def nocutconditions(self, lossper):

        # CONDIÇÃO 1 PERDAS MÉDIAS % JOULE ABAIXO DO VALOR MÁXIMO
        meanlossjoule = self.df['% JL'].mean()
        print('\nCONDITION 1')
        if meanlossjoule <= lossper:
            print('OK. THE AC LOSS MEAN {0:.2f} % IS <= {1:.2f} %.'.format(
                meanlossjoule, lossper))
        else:
            print('NG. THE AC LOSS MEAN {0:.2f} % IS >= {1:.2f} %.'.format(
                meanlossjoule, lossper))

        self.max_curr = self.df['I(A)'].max()

        # CONDIÇÃO 2 - TEMPERATURA DO CABO ABAIXO DE 94° PARA ACSR E AAAC

        def testcond3(curr):
            df = self.df

            df['x1'] = [
                1 if i > curr else 0 for i in df[('I(A)')]]
            df['y1'] = np.zeros(df.shape[0], dtype=int)
            df['y1'] = df['x1'] * (df['x1'] + df['y1'])
            if(df['y1'].max() >= 96 or (df['y1'].groupby(df.index.year).sum() >= 431).any()):
                return 1
            else:
                return 0

        def iteration(step):
            new_curr = self.max_curr
            while not testcond3(new_curr):
                new_curr -= step
            return new_curr

        cur_longdur = iteration(step=1)

        # CALCULATING THE TEMPERATURE FOR LONG DURATION CURRENT

        # Fixing the ambient temperature to a constant value
        climavars.tamb = round(self.df['Temperature'].groupby(
            self.df.index.day).max().mean())

        temp_longdur = self.temp(cur_longdur)

        # REN 191
        x = [50, 55, 60, 64, 65, 70, 75, 80, 90]
        y = [1.42, 1.33, 1.26, 1.24, 1.23, 1.19, 1.17, 1.15, 1.12]
        # POLINOMIO 4 GRAU MELHOR COMPORTAMENTO COM CURVA
        pol = np.polyfit(x, y, 4)
        k_factor = np.polyval(pol, temp_longdur)
        cur_shortdur = cur_longdur * k_factor
        temp_shortdur = self.temp(cur_shortdur)

        # CALCULATING LONG AND SHORT-TERM CURRENTS STANDARD 738- 2006
        print('\nCONDITION 2')

        print('\nLONG-TERM CURRENT: {0:.2f} A'.format(cur_longdur))
        print(
            'CABLE TEMPERATURE IN IEEE STD 738-2006: {0:.2f} °C'.format(temp_longdur))

        # CALCULATING LONG AND SHORT-TERM CURRENTS DATABASE
        print('SHORT-TERM CURRENT: {0:.2f} A'.format(cur_shortdur))
        print(
            'CABLE TEMPERATURE IN IEEE STD 738-2006: {0:.2f} °C'.format(temp_shortdur))
        if (temp_shortdur > 94):
            print('WARNING: CABLE TEMPERATURE ABOVE 94°C. NG')

        else:
            print('OK')

        # CALCULO CORRENTE DE LONGA CONDIÇÃO 3
        print('\nCONDITION 3')

        print('\nHOURS ABOVE LONG-TERM CURRENT')
        print(self.df['y1'].groupby(self.df.index.year).sum().to_string())

        print(
            '\nLONG-TERM CURRENT THAT RESPECT THE CONDITION 3: {0:.0f} A'.format(cur_longdur))
        print(
            'SHORT-TERM CURRENT (REN 191): {0:.0f}\nMAXIMUM CURRENT: {1:.0f} A'.format(cur_shortdur, self.max_curr))
        if cur_shortdur > self.max_curr:
            print('OK')
        else:
            print('NG')
        print('\f', end='')

    def curvecur(self, nameoutput):
        df = self.df  # TODO
        tpd_mean = df['I(A)'].groupby(df.index.hour).mean()
        tpd_max = df['I(A)'].groupby(df.index.hour).max()

        sns.set()
        base_color0 = sns.color_palette()[0]
        base_color1 = sns.color_palette()[1]

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax2 = ax.twinx()

        # Major ticks every 20, minor ticks every 5
        major_ticks = np.arange(0, tpd_mean.max(), 100)
        minor_ticks = np.arange(0, tpd_mean.max(), 50)

        day_major_ticks = np.arange(0, len(tpd_mean.index), 5)
        day_minor_ticks = np.arange(0, len(tpd_mean.index), 1)

        ax.set_xticks(day_major_ticks)
        ax.set_xticks(day_minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)

        ax.grid(which='major', alpha=0.5)
        ax.grid(which='minor', alpha=0.2)

        max_value = df['I(A)'].max()

        ax.set_ylim(0, max_value * (1.1))
        ax2.set_ylim(0, max_value * (1.1))

        ax.bar(tpd_mean.index, tpd_mean.values,
               width=0.8, color=base_color0)
        ax2.bar(tpd_max.index, (tpd_max.values - tpd_mean.values),
                width=0.8, bottom=tpd_mean.values, color=base_color1)

        ax.set_xlabel('Hour')
        ax.set_ylabel('I(A)')

        ax.legend(['Mean I(A)'], loc='upper left')
        ax2.legend(['Maximum I(A)'], loc='upper right')

        cwd = Path.cwd()
        imagename = '{}.png'.format(nameoutput)
        imagefilepath = cwd / imagename
        if imagefilepath.exists():
            pass
        else:
            plt.savefig('{}.png'.format(nameoutput))

    def curvecurtemp(self, nameoutput):
        df = self.df
        currmean = df['I(A)'].groupby(df.index.hour).mean()
        currmax = df['I(A)'].groupby(df.index.hour).max()

        tempmean = df['Cable Temp'].groupby(df.index.hour).mean()
        tempmax = df['Cable Temp'].groupby(df.index.hour).max()

        sns.set()
        base_color0 = sns.color_palette()[0]
        base_color1 = sns.color_palette()[1]
        base_color2 = sns.color_palette()[2]
        base_color3 = sns.color_palette()[3]

        fig = plt.figure(figsize=(10, 5))
        # CURRENT
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = ax1.twinx()
        # TEMPERATURE
        ax3 = fig.add_subplot(1, 2, 2)
        ax4 = ax3.twinx()

        max_curr = df['I(A)'].max()
        max_temp = df['Cable Temp'].max()

        ax1.set_ylim(0, max_curr * 1.1)
        ax2.set_ylim(0, max_curr * 1.1)
        ax2.set_yticks([])

        ax3.set_ylim(0, max_temp * 1.1)
        ax4.set_ylim(0, max_temp * 1.1)
        ax4.set_yticks([])

        bar1 = ax1.bar(currmean.index, currmean.values,
                       width=0.8, color=base_color0, label='Mean I(A)')
        bar2 = ax2.bar(currmax.index, (currmax.values - currmean.values),
                       width=0.8, bottom=currmean.values, color=base_color1, label='Max I(A)')

        bar3 = ax3.bar(tempmean.index, tempmean.values,
                       width=0.8, color=base_color2, label='Mean Tcond(°C)')
        bar4 = ax4.bar(tempmax.index, (tempmax.values - tempmean.values),
                       width=0.8, bottom=tempmean.values, color=base_color3, label='Max Tcond(°C)')

        # added these three lines
        bars_curr = [bar1, bar2]
        temp_curr = [bar3, bar4]

        labs1 = [l.get_label() for l in bars_curr]
        ax1.legend(bars_curr, labs1, loc='upper left')

        labs2 = [l.get_label() for l in temp_curr]
        ax3.legend(temp_curr, labs2, loc='upper left')

        # plt.tight_layout()

        cwd = Path.cwd()
        imagename = '{}.png'.format(nameoutput)
        imagefilepath = cwd / imagename
        if imagefilepath.exists():
            pass
        else:
            plt.savefig('{}.png'.format(nameoutput))

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
