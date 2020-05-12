import numpy as np
import pandas as pd

import thermalstd
import dataclima
import solarpower


db_cable = 'DB_cables.xlsx'
csvfile = r'D:\Analise_Dados_Solares\UFV Rio do Peixe\Séries de longo prazo (Helio-Clim3)\SAO_JOAO_DO_RIO_DO_PEIXE_HC3-METEO_hour_lat-6.725_lon-38.454_2004-02-01_2019-01-30_hz1.csv'

# dictstudy_ACSR = {'Type': ['ACSR', 'ACSR', 'ACSR', 'ACSR', 'ACSR', 'ACSR'],
#                   'Name': ['Partridge', 'Linnet', 'Ibis', 'Hawk', 'Dove', 'Grosbeak']}
# cablesstudy_ACSR = pd.DataFrame(dictstudy_ACSR)

# dictstudy_AAAC = {'Type': ['AAAC_1120', 'AAAC_1120', 'AAAC_1120', 'AAAC_1120', 'AAAC_1120', 'AAAC_1120'],
#                   'Name': ['Krypton', 'Lutetium', 'Neon', 'Nitrogen', 'Nobelium', 'Oxygen']}
# cablesstudy_AAAC = pd.DataFrame(dictstudy_AAAC)


# LOADING HC3 FILE
dataclima = dataclima.helioclim3(csvfile, 'rdp.pkl')
df1 = dataclima.loading()

print('ORIGINAL DATAFRAME')
print(df1.head())

# CLIMATE VARIABLES
climavars = thermalstd.climavars(vw=1.0,
                                 em=0.5,
                                 ab=0.5,
                                 tamb=40,
                                 zl=90,
                                 lat=-6,
                                 atm=1,
                                 he=100,
                                 phi=90,
                                 hour=11,
                                 nday=172)


def study_cables(dictcables, df1):

    for i in range(dictcables.shape[0]):
        cable_type = dictcables.iloc[i, 0]
        cable_name = dictcables.iloc[i, 1]
        cablevars = thermalstd.cablevars(db_cable=db_cable,
                                         cable_type=cable_type,
                                         cable_name=cable_name)
        calc = thermalstd.Std7382006(climavars=climavars, cablevars=cablevars)
        calc.graphcable()

        # CALCULATING GROSS AND NET PRODUTION

        # Fatores de perdas considerados no cálculo da Produção de Energia (%)
    # print(df1.head())

        dataloss = solarpower.energycalc(df=df1,
                                         horizon=0.2,
                                         shadings=1.9,
                                         iam=1.4,
                                         soiling=1.5,
                                         lowirradeff=0.3,
                                         temperatureloss=10.1,
                                         modulequality=0.2,
                                         lid=2.1,
                                         mismatch=0.6,
                                         ohmicdcloss=1.1,
                                         inverterloss=1.4,
                                         plantcontroller=2.5,
                                         transf_lv_mv=1.2,
                                         transf_mv_hv=0.6,
                                         auxloadsloss=0.3,
                                         ohmicac_poi=1.3,
                                         systemunavailability=0.8,
                                         gridunavailability=0.2)

        # print(df1.head())

        '''
            Características UFV
            Modelo módulo: TSM-370DE14A(II) (380W)
            Dimensão módulo: 1960 × 992 × 40 mm
            https://www.civicsolar.com/question/how-do-you-calculate-solar-panel-efficiency

            modulearea = 1.96 * 0.992  # m²
            '''

        dfproduction = dataloss.production(modulearea=1.94432,
                                           totalpower=80.256e6,
                                           modulepower=380,
                                           trackeradd=1.276)

        # CUTTING MAX PRODUCTION IN SOLAR POWER PLANT

        linevars = thermalstd.linevars(dfproduction=dfproduction,
                                       voltage=69,
                                       powerfactor=0.95,
                                       climavars=climavars,
                                       cablevars=cablevars,
                                       extline=14,
                                       maxnetprod=61)

        # WITHOUT CUTTING MAX PRODUCTION IN SOLAR POWER PLANT

        # linevars = thermalstd.linevars(dfproduction=dfproduction,
        #                                voltage=69,
        #                                powerfactor=0.95,
        #                                climavars=climavars,
        #                                cablevars=cablevars,
        #                                extline=14)

        # outnetlimit = power value at the destiny
        # maxnetprod = max power valer at the origin

        # CONDITIONS

        dataanalysis = thermalstd.analysis(climavars=climavars,
                                           cablevars=cablevars,
                                           linevars=linevars,
                                           savexlsx=True)

        print('PLOTTING MEAN AND MAX CURRENT BARS')
        dataanalysis.curvecur('Current_Bars')

        print('PLOTTING CURRENT X TEMP BARS')
        dataanalysis.curvecurtemp('Current_Temp_Bars')

        dataanalysis.conditions()


# study_cables(cablesstudy_ACSR, df1)
# study_cables(cablesstudy_AAAC, df1)

dictstudy = {'Type': ['AAAC_1120'],
             'Name': ['Nitrogen']}
cablesstudy = pd.DataFrame(dictstudy)

study_cables(cablesstudy, df1)
