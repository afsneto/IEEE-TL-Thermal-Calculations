import numpy as np
import pandas as pd

import thermalstd as tstd
from dataclima import helioclim3
from solarpower import energycalc
# from cablepower import losspower

vw = 0.6
em = 0.5
ab = 0.5
tamb = 40
zl = 90
lat = 30
atm = 1
he = 100
phi = 90
hour = 11
nday = 172

db_cable = 'DB_cables.xlsx'
cable_type = 'ACSR'
cable_name = 'Drake'

csvfile = r'D:\Analise_Dados_Solares\UFV Rio do Peixe\Séries de longo prazo (Helio-Clim3)\SAO_JOAO_DO_RIO_DO_PEIXE_HC3-METEO_hour_lat-6.725_lon-38.454_2004-02-01_2019-01-30_hz1.csv'
voltage = 69e3
extline = 14  # km

climavars = tstd.climavars(vw=vw,
                           em=em,
                           ab=ab,
                           tamb=tamb,
                           zl=zl,
                           lat=lat,
                           atm=atm,
                           he=he,
                           phi=phi,
                           hour=hour,
                           nday=nday)

cablevars = tstd.cablevars(db_cable=db_cable,
                           cable_type=cable_type,
                           cable_name=cable_name)

calc = tstd.Std7382006(climavars=climavars, cablevars=cablevars)

# tc = 100
# ivalue = calc.current(tc)
# calc.graphcable()


# LOADING HC3 FILE

dataclima = helioclim3(csvfile, 'rdp.pkl')
df = dataclima.dfloaded()

# CALCULATING GROSS AND NET PRODUTION

'''
Fatores de perdas considerados no cálculo da Produção de Energia (%)
'''
horizon = 0.2
shadings = 1.9
iam = 1.4
soiling = 1.5
lowirradeff = 0.3
temperatureloss = 10.1
modulequality = 0.2
lid = 2.1
mismatch = 0.6
ohmicdcloss = 1.1
inverterloss = 1.4
plantcontroller = 2.5
transf_lv_mv = 1.2
transf_mv_hv = 0.6
auxloadsloss = 0.3
ohmicac_poi = 1.3
systemunavailability = 0.8
gridunavailability = 0.2


'''
Características UFV
Modelo módulo: TSM-370DE14A(II) (370W)
Dimensão módulo: 1960 × 992 × 40 mm
https://www.civicsolar.com/question/how-do-you-calculate-solar-panel-efficiency
'''
modulearea = 1.96 * 0.992  # m²
totalpower = 80.256e6
modulepower = 380  # Wp
trackeradd = 1.276

dataloss = energycalc(df=df,
                      horizon=horizon,
                      shadings=shadings,
                      iam=iam,
                      soiling=soiling,
                      lowirradeff=lowirradeff,
                      temperatureloss=temperatureloss,
                      modulequality=modulequality,
                      lid=lid,
                      mismatch=mismatch,
                      ohmicdcloss=ohmicdcloss,
                      inverterloss=inverterloss,
                      plantcontroller=plantcontroller,
                      transf_lv_mv=transf_lv_mv,
                      transf_mv_hv=transf_mv_hv,
                      auxloadsloss=auxloadsloss,
                      ohmicac_poi=ohmicac_poi,
                      systemunavailability=systemunavailability,
                      gridunavailability=gridunavailability)

dfproduction = dataloss.production(modulearea=modulearea,
                                   totalpower=totalpower,
                                   modulepower=modulepower,
                                   trackeradd=trackeradd)

# print(dfproduction.head())

# CALCULATING ELECTRIC LIMITS

voltage = 69e3
powerfactor = 0.95


dataelectric = tstd.losspower(dfproduction=dfproduction,
                              voltage=voltage,
                              powerfactor=powerfactor,
                              climavars=climavars,
                              cablevars=cablevars)

print(dataelectric.completedf().head(24))
