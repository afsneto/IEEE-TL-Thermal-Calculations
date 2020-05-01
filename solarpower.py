from calendar import month_name, monthrange
from pathlib import Path, PureWindowsPath


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import math as mt

from dataclima import helioclim3


class energycalc:
    def __init__(self, df,
                 horizon,
                 shadings,
                 iam,
                 soiling,
                 lowirradeff,
                 temperatureloss,
                 modulequality,
                 lid,
                 mismatch,
                 ohmicdcloss,
                 inverterloss,
                 plantcontroller,
                 transf_lv_mv,
                 transf_mv_hv,
                 auxloadsloss,
                 ohmicac_poi,
                 systemunavailability,
                 gridunavailability):

        self.horizon = horizon
        self.shadings = shadings
        self.iam = iam
        self.soiling = soiling
        self.lowirradeff = lowirradeff
        self.temperatureloss = temperatureloss
        self.modulequality = modulequality
        self.modulequality = modulequality
        self.lid = lid
        self.mismatch = mismatch
        self.ohmicdcloss = ohmicdcloss
        self.inverterloss = inverterloss
        self.plantcontroller = plantcontroller
        self.transf_lv_mv = transf_lv_mv
        self.transf_mv_hv = transf_mv_hv
        self.auxloadsloss = auxloadsloss
        self.ohmicac_poi = ohmicac_poi
        self.systemunavailability = systemunavailability
        self.gridunavailability = gridunavailability

        self.df = self._datatreat(df)

    def _datatreat(self, df):
        # df.drop(columns=['Top of Atmosphere', 'Code', 'Relative Humidity',
        #                  'Wind direction', 'Rainfall', 'Snowfall',
        #                  'Snow depth'],
        #         inplace=True)
        df.drop(columns=['Top of Atmosphere',
                         'Code', 'Relative Humidity',
                         'Wind direction', 'Rainfall', 'Snowfall',
                         'Snow depth'])

        return df

    def perfratio(self):

        lossesModuleFactors = {'Horizon': self.horizon,
                               'Shadings': self.shadings,
                               'IAM': self.iam,
                               'Soiling': self.soiling}

        lossesLocalFactors = {'Low Irradiance efficiency fall-off': self.lowirradeff,
                              'Temperature': self.temperatureloss,
                              'Module Quality': self.modulequality,
                              'LID': self.lid,
                              'Mismatch': self.mismatch,
                              'Ohmic (DC)': self.ohmicdcloss,
                              'Inverter': self.inverterloss,
                              'Plant Controller': self.plantcontroller,
                              'Transformers LV-MV': self.transf_lv_mv,
                              'Transformers MV-HV': self.transf_mv_hv,
                              'Auxiliary Loads': self.auxloadsloss,
                              'Ohmic AC until POI': self.ohmicac_poi,
                              'System Unavailability': self.systemunavailability,
                              'Grid Unavailability': self.gridunavailability}

        '''
        Performace Ratio (Desempenho global)
        '''
        dflocalloss = pd.DataFrame(data=lossesLocalFactors.values(),
                                   index=lossesLocalFactors.keys(),
                                   columns=['value'])

        dfmodloss = pd.DataFrame(data=lossesModuleFactors.values(),
                                 index=lossesModuleFactors.keys(),
                                 columns=['value'])

        vector1 = np.array(dflocalloss['value'])
        frac1 = (1 - (vector1)/100)
        local_losses = np.cumprod(frac1, dtype=float)[-1]

        vector2 = np.array(dfmodloss['value'])
        frac2 = (1 - (vector2)/100)
        module_losses = np.cumprod(frac2, dtype=float)[-1]

        perf_ratio = local_losses * module_losses
        # print('Desempenho global: {0:.2f} %'.format(perf_ratio * 100))
        return module_losses, local_losses

    def production(self, modulearea, totalpower, modulepower, trackeradd):

        module_losses, local_losses = self.perfratio()
        # Eficiencia do painel solar
        eff = modulepower / (modulearea * 1000)

        # Calculo total area
        number_modul = totalpower / modulepower  # MWp total / MWp
        totalmodulearea = number_modul * modulearea  # m²

        # Irradiância da UFV no plano dos módulos (kWh/m²)

        irrad_eff = self.df['Global Horiz'] * \
            trackeradd * module_losses

        # Calculating gross and net production

        gross_production = irrad_eff * totalmodulearea * eff / 1e6  # Wh para MWh
        # net_prodution = gross_production * perf_ratio
        net_prodution = gross_production * local_losses
        dfprod = pd.DataFrame({'Gross Prod': gross_production,
                               'Net Prod': net_prodution})

        dfconcat = pd.concat([self.df, dfprod], axis=1, sort=False)
        # https://pandas-docs.github.io/pandas-docs-travis/user_guide/merging.html

        return dfconcat
