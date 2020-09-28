import dask.dataframe as dd
import numpy as np
from scipy.signal import welch
import matplotlib.pyplot as plt


def att():
    """ This function is designed to plot the images about attitude
    """
    dd_att = dd.read_csv(urlpath='..//output//taiji-01-0811-attitude-2019-09.txt', sep='\t', header=None,
                         engine='c', skiprows=44, storage_options=dict(auto_mkdir=False),
                         names=['gps_time', 'gcrs_x', 'gcrs_y', 'gcrs_z', 'tf_x', 'tf_y', 'tf_z'],
                         dtype=np.float64, encoding='gb2312')
    
    freq_igrf_x, psd_igrf_x = welch(
        dd_att.gcrs_x.compute().to_numpy(), 1., 'hanning', dd_att.__len__(), scaling='density')
    freq_igrf_y, psd_igrf_y = welch(
        dd_att.gcrs_y.compute().to_numpy(), 1., 'hanning', dd_att.__len__(), scaling='density')
    freq_igrf_z, psd_igrf_z = welch(
        dd_att.gcrs_z.compute().to_numpy(), 1., 'hanning', dd_att.__len__(), scaling='density')

    plt.style.use(['science', 'no-latex', 'high-vis'])
    fig, ax = plt.subplots(figsize=(15, 8))
    plt.plot(dd_att.gps_time.compute().to_numpy(), dd_att.gcrs_x.compute().to_numpy(),
             linewidth=2, label='gcrs2srf_x')
    plt.plot(dd_att.gps_time.compute().to_numpy(), dd_att.gcrs_y.compute().to_numpy(),
             linewidth=2, label='gcrs2srf_y')
    plt.plot(dd_att.gps_time.compute().to_numpy(), dd_att.gcrs_z.compute().to_numpy(),
             linewidth=2, label='gcrs2srf_z')
    plt.tick_params(labelsize=25, width=2.9)
    ax.yaxis.get_offset_text().set_fontsize(24)
    ax.xaxis.get_offset_text().set_fontsize(24)
    plt.xlabel('GPS time [s]', fontsize=20)
    plt.ylabel('Attitude angle [degree]', fontsize=20)
    plt.legend(fontsize=20, loc='best')
    plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    
    fig, ax = plt.subplots(figsize=(15, 8))
    plt.loglog(freq_igrf_x, np.sqrt(psd_igrf_x),
               linewidth=2, label='gcrs2srf_x')
    plt.loglog(freq_igrf_y, np.sqrt(psd_igrf_y),
               linewidth=2, label='gcrs2srf_y')
    plt.loglog(freq_igrf_z, np.sqrt(psd_igrf_z),
               linewidth=2, label='gcrs2srf_z')
    plt.tick_params(labelsize=25, width=2.9)
    ax.yaxis.get_offset_text().set_fontsize(24)
    ax.xaxis.get_offset_text().set_fontsize(24)
    plt.xlabel('Frequency [Hz]', fontsize=20)
    plt.ylabel(r'$PSD \quad [deg/\sqrt{hz}]$', fontsize=20)
    plt.legend(fontsize=20, loc='best')
    plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    
    plt.show()
    
    
def atmos():
    dd_atmos = dd.read_csv(urlpath='..//output//taiji-01-0860-air-density-2019-09.txt', header=None,
                           engine='c', skiprows=1, storage_options=dict(auto_mkdir=False),
                           names=['air_density'], dtype=np.float64, encoding='gb2312')
    
    plt.style.use(['science', 'no-latex', 'high-vis'])
    fig, ax = plt.subplots(figsize=(15, 8))
    plt.plot(dd_atmos.air_density.compute().to_numpy(), linewidth=2)
    plt.tick_params(labelsize=25, width=2.9)
    ax.yaxis.get_offset_text().set_fontsize(24)
    ax.xaxis.get_offset_text().set_fontsize(24)
    plt.xlabel('GPS time [s]', fontsize=20)
    plt.ylabel(r'$Air \quad density [kg/m^3]$', fontsize=20)
    plt.legend(fontsize=20, loc='best')
    plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    
    plt.show()
    
    
def gcrs():
    dd_gcrs = dd.read_csv(urlpath='..//output//taiji-01-0866-gcrs-2019-09.txt', header=None,
                          engine='c', skiprows=50, storage_options=dict(auto_mkdir=False), sep='\t',
                          names=['gps_time', 'xpos', 'ypos', 'zpos', 'xrvel', 'yrvel', 'zrvel'],
                          dtype=np.float64, encoding='gb2312')
    
    plt.style.use(['science', 'no-latex', 'high-vis'])
    fig, ax = plt.subplots(figsize=(15, 8))
    plt.plot(dd_gcrs.xpos.compute().to_numpy(), linewidth=2, label='xpos_gcrs')
    plt.plot(dd_gcrs.ypos.compute().to_numpy(), linewidth=2, label='ypos_gcrs')
    plt.plot(dd_gcrs.zpos.compute().to_numpy(), linewidth=2, label='zpos_gcrs')
    plt.tick_params(labelsize=25, width=2.9)
    ax.yaxis.get_offset_text().set_fontsize(24)
    ax.xaxis.get_offset_text().set_fontsize(24)
    plt.xlabel('GPS time [s]', fontsize=20)
    plt.ylabel(r'$Position [km]$', fontsize=20)
    # plt.legend(fontsize=20, loc='best')
    plt.legend(fontsize=20, loc='lower left', bbox_to_anchor=(0, 1, 1, .1), ncol=3, mode='expand')
    plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    
    fig, ax = plt.subplots(figsize=(15, 8))
    plt.plot(dd_gcrs.xrvel.compute().to_numpy(), linewidth=2, label='xrvel_gcrs')
    plt.plot(dd_gcrs.yrvel.compute().to_numpy(), linewidth=2, label='yrvel_gcrs')
    plt.plot(dd_gcrs.zrvel.compute().to_numpy(), linewidth=2, label='zrvel_gcrs')
    plt.tick_params(labelsize=25, width=2.9)
    ax.yaxis.get_offset_text().set_fontsize(24)
    ax.xaxis.get_offset_text().set_fontsize(24)
    plt.xlabel('GPS time [s]', fontsize=20)
    plt.ylabel(r'$Relative \quad velocities [km/s]$', fontsize=20)
    # plt.legend(fontsize=20, loc='best')
    plt.legend(fontsize=20, loc='lower left', bbox_to_anchor=(
        0, 1, 1, .1), ncol=3, mode='expand')
    plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)

    plt.show()
    
    
if __name__ == '__main__':
    # att()
    # atmos()
    gcrs()
