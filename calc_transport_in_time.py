'''
Find transport to Port Aransas region in time from previously-existing files,
calculated in find_coastal_path_connectivity.py.
'''

from glob import glob
import numpy as np
import pandas as pd

# from make_conn_plots
# names = ['bpm', 'bpa', 'bpoc', 'bgalv', 'batch', 'bterr', 'bbara', 'all']
# # distances = [400, 555, 645, 840, 1175, 1280, 1350]  # atch 1190, bterr1290
# boxes = [np.arange(76,88), np.arange(107,119), np.arange(126,138), np.arange(165,177),
#          np.arange(234,246), np.arange(257,269), np.arange(271, 283), np.arange(0,342)]
boxes = {'bpm': np.arange(76,88), 'bpa': np.arange(107,119),
         'bpoc': np.arange(126,138), 'bgalv': np.arange(165,177),
         'batch': np.arange(234,246), 'bterr': np.arange(257,269),
         'bbara': np.arange(271, 283), 'all': np.arange(0,342) }

def calc_df(from='all', to='bpa'):
    '''Calculate dataframes of transport from 'from' to 'to', in time.

    Make a column for each simulation.'''

    # # taken from plot_conn2coast_2boxes.py for Port Aransas region
    # boxes = np.arange(103,123)


    Files = np.sort(glob('calcs/alongcoastconn/20??-??.npz'))
    Files = np.sort(Files)

    # create dataframe to stick calculations into
    dstart = Files[0].split('/')[-1].split('.')[0]
    dend = Files[-1].split('/')[-1].split('.')[0]
    dfdates = pd.date_range(start=dstart, end=dend, freq='2880S')  # 48 min (drifter output freq)
    df = pd.DataFrame(index=dfdates)


    for File in Files:
        d = np.load(File)
        inbox = d['inbox']  # coastal boxes x drifters that make it inside x box boundary crosses timing (up to 5)
        # use this to relate to the drifters overall since we are only examining a subset here
        iinside = d['iinside']
        d.close()

        # only using first entrance to box of each drifter right now. Worried that
        # if I use subsequent entrances as well that the normalization between
        # runs will lose meaning.
        inbox = inbox[boxes{to}, :, 0]

        # grab all times for when a drifter is entering these boxes.
        # don't care which boxes, from where they come; only when they enter
        # this includes multiple entries from the same drifter potentially
        i, j = np.where(~np.isnan(inbox))
        times = inbox[i, j]

        # Get signal: count drifter entering in time
        # uses number of bins for 48 min drifter time res; 0.8 is 48 min drifter output
        ndrifters, bins = np.histogram(times, bins=int(30*(24/.8)))

        # # normalize by how many drifters reached the coastline for this simulation,
        # # so results are percent drifters near coastline reaching Port Aransas region
        # inbox /= iinside.sum()


        dstart = File.split('/')[-1].split('.')[0]
        # dend = File.split('/')[-1].split('.')[0]
        dfdates = pd.date_range(start=dstart, periods=bins.size-1, freq='2880S')  # 48 min (drifter output freq)
        # dftemp holds number drifters entering boxes in time for simulation
        dftemp = pd.DataFrame(index=dfdates, data={dstart: ndrifters})
            # dftemp = pd.DataFrame(index=netCDF.num2date(tp, 'seconds since 1970-01-01  00:00:00'), data={simstartdate: numexit})
        # need to normalize across columns at the end since different number of drifters
        # available to be entering boxes at different times
        df = df.join(dftemp)  # add column to dataframe

        df.to_csv('calcs/alongcoastconn/df_2013-0508.csv')  # save every time step


def combine_cols():
    '''Combine columns of simulations from calc_df for easy analysis.'''

    # load in df
    df = pd.read_csv('calcs/alongcoastconn/df_2013-0508.csv', parse_dates=True, index_col=0)

    nfiles = (~np.isnan(df)).astype(int).sum(axis='columns')
    y = df.sum(axis='columns').divide(nfiles)
    # y = df.sum(axis='columns').divide(nfiles).resample('15min', base=0).interpolate()

    dfnew = pd.DataFrame(index=df.index, data={'drifters': y})
    dfnew.to_csv('calcs/alongcoastconn/df_2013-0508_combined.csv')
