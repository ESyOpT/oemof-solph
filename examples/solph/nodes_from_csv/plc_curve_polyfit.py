# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import scipy.stats as sp
import matplotlib
import matplotlib.pyplot as plt

# global plotting options
plt.rcParams.update(plt.rcParamsDefault)
matplotlib.style.use('ggplot')
plt.rcParams['lines.linewidth'] = 1.2
plt.rcParams['axes.facecolor'] = 'silver'
plt.rcParams['xtick.color'] = 'k'
plt.rcParams['ytick.color'] = 'k'
plt.rcParams['text.color'] = 'k'
plt.rcParams['axes.labelcolor'] = 'k'
plt.rcParams.update({'font.size': 10})
plt.rcParams['image.cmap'] = 'Spectral'

# read file
file = ('results/'
        'results_dispatch_prices_DE_2016-07-25 15:30:49.914634_nep_2014_aggr'
        '.csv')

df_raw = pd.read_csv(file, parse_dates=[0], index_col=0, keep_date_col=True)
df_raw.head()
df_raw.columns


# %% polynom fitting

# prepare dataframe for fit
residual_load = df_raw['DE_load'] + df_raw['AT_load'] + df_raw['LU_load'] - \
                df_raw['DE_wind'] - df_raw['AT_wind'] - df_raw['LU_wind'] - \
                df_raw['DE_solar'] - df_raw['AT_solar'] - df_raw['LU_solar']
df = pd.concat([residual_load, df_raw['eex_day_ahead_2014'],
                df_raw['power_price_model']], axis=1)
df.columns = ['res_load', 'price_real', 'price_model']

# fit polynom of 3rd degree
z = np.polyfit(df['res_load'], df['price_real'], 3)
p = np.poly1d(z)

# save and plot results
df['price_polynom'] = p(df['res_load'])

df['residuals'] = df['price_real'] - \
                  df['price_model']

# %% create normal distributed volatility

# param[0] and param[1] are the mean and
# the standard deviation of the fitted distribution
mu, sigma = sp.norm.fit(df['residuals'])

df['random_norm'] = np.random.normal(mu, sigma, 8760)

df['price_model_volatility'] = df['price_model'] + \
                               df['random_norm']


# %% spread analysis

# data
df_spread = pd.DataFrame()

df_spread['price_real'] = df['price_real']

df_spread['spread_3h'] = df['price_real'].resample('3h').max() - \
    df['price_real'].resample('3h').min()

df_spread['spread_6h'] = df['price_real'].resample('6h').max() - \
    df['price_real'].resample('6h').min()

df_spread['spread_12h'] = df['price_real'].resample('12h').max() - \
    df['price_real'].resample('12h').min()

df_spread['spread_24h'] = df['price_real'].resample('24h').max() - \
    df['price_real'].resample('24h').min()

df_spread['spread_48h'] = df['price_real'].resample('48h').max() - \
    df['price_real'].resample('48h').min()

df_spread['spread_96h'] = df['price_real'].resample('96h').max() - \
    df['price_real'].resample('96h').min()

df_spread['spread_192h'] = df['price_real'].resample('192h').max() - \
    df['price_real'].resample('192h').min()


# %% plotting

#df.plot(kind='scatter', x='res_load', y='price_real')

#df.plot(kind='scatter',
#        x='price_real', y='price_model')

#df[:][['price_real', 'price_model']].plot(linewidth=1.2, subplots=True,
#                                          drawstyle='steps',
#                                          color=['grey', 'r', 'b'],
#                                          ylim=[-100, 100])

#residuals = pd.DataFrame()
#residuals = pd.concat([df['residuals'][0:2190],
#                       df['residuals'][2190:4380],
#                       df['residuals'][4380:6570],
#                       df['residuals'][6570:8760]], axis=1)
#residuals.plot.hist(bins=100, subplots=True, legend=None)
#residuals.columns=['Q1', 'Q2', 'Q3', 'Q4']
#

df[['residuals', 'random_norm']].plot.hist(bins=100, subplots=True,
                                           sharex=True, sharey=True,
                                           layout=(1,2))

plt.show()

#df[0:24 * 31][['price_real', 'price_model',
#               'price_model_volatility']].plot(linewidth=1.2,
#                                               subplots=True,
#                                               drawstyle='steps',
#                                               ylim=[-100, 100])
#
#plt.show()

#fig, axes = plt.subplots(nrows=7, sharey=True)
#fig.suptitle('Spread nach Zeitintervall', fontsize=16)
#
#df_spread[['spread_3h']].dropna().plot(kind='line', drawstyle='steps', ax=axes[0])
#df_spread[['spread_6h']].dropna().plot(kind='line', drawstyle='steps', ax=axes[1])
#df_spread[['spread_12h']].dropna().plot(kind='line', drawstyle='steps', ax=axes[2])
#df_spread[['spread_24h']].dropna().plot(kind='line', drawstyle='steps', ax=axes[3])
#df_spread[['spread_48h']].dropna().plot(kind='line', drawstyle='steps', ax=axes[4])
#df_spread[['spread_96h']].dropna().plot(kind='line', drawstyle='steps', ax=axes[5])
#df_spread[['spread_192h']].dropna().plot(kind='line', drawstyle='steps', ax=axes[6])
#
#plt.show()