from scipy.signal import kaiserord, firwin, filtfilt, lfilter, freqz, welch

import numpy as np
import dask.dataframe as dd

import matplotlib.pyplot as plt


def main():
    dd_gps = dd.read_csv(urlpath='..//output//09//taiji-01-0860-earth-fixed-system-2019-09-after-fortran.txt', sep='\s+', header=None,
                         engine='c', skiprows=91, storage_options=dict(auto_mkdir=False),
                         names=['gps_time', 'rcv_time', 'pos_e1', 'pos_e2', 'pos_e3',
                                'vel1', 'vel2', 'vel3', 'acc1', 'acc2', 'acc3',
                                'vac_qualflg', 'vel_diffed1', 'vel_diff2', 'vel_diffed3',
                                'pj_qualflg', 'vp_qualflg'],
                         dtype=np.float64, encoding='gb2312')

    # the desired width of the transition from pass to stop
    width = 0.05 / 0.5
    # the desired attenuation in the stop band in db: ripple_db
    # compute the kaiser parameter for the fir filter
    n, beta = kaiserord(80.0, width)
    print('The length of the lowpass filter is', n, '.')
    # use firwin with a kaiser window
    taps = firwin(n, 0.01, window=('kaiser', beta),
                  pass_zero='lowpass', nyq=0.5)
    # use filtfilt to filter x with the fir filter
    flags = ['pos_e1', 'pos_e2', 'pos_e3', 'vel1', 'vel2', 'vel3']
    filtered_list= np.zeros([dd_gps.__len__(), flags.__len__()])
    for index, flag in enumerate(flags):
        filtered_list[:, index] = filtfilt(taps, 1.0, dd_gps[flag].compute().to_numpy())

    np.savetxt(
        '..//output//09//taiji-01-0860-earth-fixed-system-2019-09-after-fortran-filtered.txt',
        filtered_list, fmt='%.18e', newline='\n', header='xpos ypos zpos xvel yvel zvel')


if __name__ == '__main__':
    main()
