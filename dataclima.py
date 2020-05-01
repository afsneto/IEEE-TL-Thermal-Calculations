import numpy as np
import pandas as pd
from pandas.tseries import converter
from pathlib import Path

from tqdm import tqdm
from datetime import datetime
from calendar import monthrange
from calendar import month_name

import matplotlib.pyplot as plt
import seaborn as sns

import swifter
import calendar
import pytz


class helioclim3:
    def __init__(self, hc3csvfile, pkl_filename):
        self.csvfile = hc3csvfile  # csv file

        cwd = Path.cwd()
        pkl_file = cwd / pkl_filename

        def savepandasfile(df, output=pkl_filename):
            df = self.loading()
            df.to_pickle(output)

        def loadpandasfile(file=pkl_filename):
            try:
                return pd.read_pickle(file)
            except:
                print('LOADING FILE ERROR\n')

        if pkl_file.exists():
            print('LOADING .PKL FILE\n')
            self.df = loadpandasfile()
        else:
            savepandasfile(self.csvfile, pkl_filename)
            self.df = loadpandasfile()

    def dfloaded(self):
        return self.df

    def _fixhours(self, str_datetime):
        # fix hour format 24:00:00
        if '24:' in str_datetime:
            str_datetime = str_datetime.replace(
                ':00', '', 1).replace('24:', '00:')
            return pd.to_datetime(str_datetime, format='%d/%m/%Y %H:%M')
        else:
            return pd.to_datetime(str_datetime, format='%d/%m/%Y %H:%M')

    def loading(self, datafix=True):
        df = pd.read_csv(self.csvfile, skiprows=31, sep=";",
                         parse_dates=[['# Date', 'Time']])

        # Rename columns
        df.rename(columns={"# Date_Time": "Date"}, inplace=True)

        # Temperature °F to °C
        df["Temperature"] = [(i - 273.15)
                             for i in df["Temperature"]]
        # Fix date
        # tqdm.pandas()
        df.loc[:, "Date"] = df.Date.swifter.apply(self._fixhours)

        # Assign "Date" column to index
        df.set_index("Date", inplace=True)
        # Convert Date to correct timezone
        # https://en.wikipedia.org/wiki/List_of_tz_database_time_zones#List
        # https://stackoverflow.com/questions/22800079/converting-time-zone-pandas-dataframe
        buenos_aires = pytz.timezone('America/Buenos_Aires')
        df.index = df.index.tz_localize(pytz.utc).tz_convert(buenos_aires)

        # Fix values -999
        if datafix:
            for i in df.columns:
                df.loc[df[str(i)] <= -999, str(i)] = 0
        else:
            pass

        return df

    def plot(self, date1: str, date2: str, data='Global Horiz'):
        df = self.loading()
        converter.register()

        d1 = datetime.strptime(date1, '%Y-%m-%d')
        d2 = datetime.strptime(date2, '%Y-%m-%d')

        df2plot = df.loc[d1:d2]
        sns.set(style="darkgrid")
        f, ax = plt.subplots(figsize=(10, 5))
        sns.lineplot(x=df2plot.index, y=df2plot[data])

        # Removing top and right borders
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Finalize the plot
        sns.despine(bottom=True)
        plt.setp(f.axes, xticks=[],
                 xlabel='Interval\nfrom: {0}\nto: {1}'.format(date1, date2),
                 ylabel='Solar Irradiation (W/m²)')
        plt.tight_layout(h_pad=2)

        plt.savefig("output.png")
        print('Image saved.')

    def analysis(self):
        """Characteristics of the data and NULL values
        """
        df = self.loading(datafix=False)
        noIrrad = df[df['Global Horiz'] == -999]
        noIrradDays = noIrrad.index.date
        noDataIrrad = len(noIrradDays)
        totalIrrad = len(df['Global Horiz'])
        percDataIrrad = (noDataIrrad/totalIrrad) * 100

        yearsIrrad = sorted(set(df.index.year.values))

        print('\nIntervalo de dados de medição: {0:d} a {1:d}'.format(
            min(yearsIrrad), max(yearsIrrad)))

        print('Número de linhas sem dados de irradiação: {0}'.format(
            noDataIrrad))
        print('Número total de linhas: {0}'.format(totalIrrad))
        print('Porcentagem de linhas sem dados de irradiação: {0:2.4f} %'.format(
            percDataIrrad))

        print('\nDias do ano sem registro de irradiação:')
        for i in sorted(set(noIrradDays)):
            print(i.strftime('%d/%m/%Y'))

        code = [0, 1, 2, 5, 6]
        numberbyCode = {i: len(df[df["Code"] == i]) for i in code}
        idbyCode = {0: 'no data', 1: 'sun below horizon',
                    2: 'satellite assessment', 5: 'interpolation in time', 6: 'forecast'}

        for i in numberbyCode.keys():
            print("{0}: {1} - {2:2.1f}%".format(
                idbyCode[i], numberbyCode[i], (numberbyCode[i] / totalIrrad)*100))

        df.info().to_string()

    def averSolarIrrad(self):
        """Calculates the average values for the irradiation kW/m²
        """
        # https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#anchored-offsets

        df = self.df
        df.drop(columns=['Top of Atmosphere', 'Code', 'Relative Humidity',
                         'Wind direction', 'Rainfall', 'Snowfall',
                         'Snow depth', 'month', 'year'], inplace=True)

        print(df.head().to_string())

        # Calculo da media para cada mes
        Wh_m2_mes = df.groupby(df.index.month).mean()['Global Horiz']
        kWh_m2_dia = df.groupby(df.index.month).mean()[
            'Global Horiz'] * 24/1000
        Months = Wh_m2_mes.index.to_list()

        result = {'kWh/m²_Diário': kWh_m2_dia,
                  'Wh/m²_Mensal': Wh_m2_mes, 'Month': Months}
        dfIrrad = pd.DataFrame(result)
        print(dfIrrad)
        print('Média diária Irradiação Solar: {0:.2f} kWh/m²/dia'.format(
            dfIrrad.loc[:, 'kWh/m²_Diário'].mean()))

    def savecsv_mean_month_yearly(self, zerovalues=True):
        """Saves a .csv file with the mean mensal of each year (kWh/m²_Month)
        and a pd.describe() of this .csv.

        Keyword Arguments:
            zerovalues {bool} -- [Consider or not the zero values of solar
            irradiation] (default: {True})
        """
        if zerovalues:
            df = self.df
        else:
            df = self.df
            df = df.loc[df['Global Horiz'] != 0]

        month_names = [month_name[i] for i in range(1, 13)]
        df_all_months = pd.DataFrame(columns=month_names)

        for i in range(1, 13):
            df_month_by_year = self.filter_df_month(i)
            df_all_months[month_name[i]] = df_month_by_year['Global Horiz']

        df_all_months.describe().to_csv('Calculos_Media_Mensal.csv')
        df_all_months.to_csv('Media_Mensal.csv')

    def boxplot_mean_month_yearly(self, zerovalues=True):
        """Boxplot monthly in a year interval

        Keyword Arguments:
            zerovalues {bool} -- [description] (default: {True})
        """
        if zerovalues:
            df = self.df
        else:
            df = self.df
            df = df.loc[df['Global Horiz'] != 0]

        month_names = [month_name[i] for i in range(1, 13)]
        df_all_months = pd.DataFrame(columns=month_names)

        for i in range(1, 13):
            df_month_by_year = self.filter_df_month(i)
            df_all_months[month_name[i]] = df_month_by_year['Global Horiz']

        sns.set()
        plt.title('{0} - {1}'.format(min(df_all_months.index.values),
                                     max(df_all_months.index.values)))

        sns.boxplot(data=df_all_months,
                    palette='Blues')
        plt.xticks(rotation=90)

        # plt.show()
        years_interval = df_all_months.index.values
        plt.savefig('boxplot{0} - {1}.png'.format(min(years_interval),
                                                  max(years_interval)), format='png')
        plt.clf()

    def save_barplot_yearly_sns(self, zerovalues=True):
        """Plot images by month along the data years

        Keyword Arguments:
            zerovalues {bool} -- [Consider or not the zero values of solar
            irradiation] (default: {True})
        """
        if zerovalues:
            df = self.df
        else:
            df = self.df
            df = df.loc[df['Global Horiz'] != 0]

        sns.set()
        for i in range(1, 13):
            df_month_by_year = self.filter_df_month(i)
            plt.title('{}'.format(calendar.month_name[i]))

            sns.barplot(x=df_month_by_year.index,
                        y='Global Horiz',
                        data=df_month_by_year,
                        color='blue',
                        capsize=.2)

            plt.xticks(rotation=90)
            plt.savefig('{}_{}.png'.format(
                i, calendar.month_name[i]), format='png')

            plt.clf()

    def save_boxplot_monthly(self, zerovalues=True):
        if zerovalues:
            df = self.df
        else:
            df = self.df
            df = df.loc[df['Global Horiz'] != 0]
        sns.set()
        plt.title('{0} - {1}'.format(min(df.index.year), max(df.index.year)))
        # month_names = [calendar.month_name[i] for i in range(1, 13)]
        # plt.xticks(np.arange(1, 13), month_names, rotation=20)
        sns.boxplot(x='month',
                    y='Global Horiz',
                    data=df,
                    palette='Blues')
        plt.savefig('{0} - {1}.png'.format(min(df.index.year),
                                           max(df.index.year)), format='png')
        plt.clf()

    def filter_df_month(self, month=1):
        df_month = self.df_monthly(month)
        df_month_by_year = df_month.groupby(df_month.index.year).mean()
        return df_month_by_year

    def df_monthly(self, month=1):
        df = self.df
        df_month = df[df['month'] == month]
        return df_month

    def describe_ghi(self):
        df = self.df
        df_by_month = df.groupby(df.index.month)
        df_ghi_month = df_by_month['Global Horiz']
        df_ghi_month.describe().to_csv('describe_ghi_month.csv')
        print(df_ghi_month.describe())
