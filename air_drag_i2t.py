import dask.dataframe as dd
import numpy as np
import matplotlib.pyplot as plt

from astropy import coordinates as coor
from astropy import units as u
from astropy.time import Time as atime
from auromat.coordinates import transform


def air_density(years, months, date_i):
    dd_aird_i = dd.read_csv(urlpath='..//output//taiji-01-0222-air-drag-gcrs-2019-09.txt', header=None,
                            engine='c', skiprows=35, storage_options=dict(auto_mkdir=False), sep='\s+',
                            names=['gps_time', 'xacc', 'yacc', 'zacc', 'xacc_s', 'yacc_s', 'zacc_s'],
                            dtype=np.float64, encoding='gb2312')
    utc_time = dd_aird_i.gps_time.compute().to_numpy()

    dates = np.floor((utc_time - utc_time[0]) / 86400.).astype(np.int) + date_i
    hours = np.floor((utc_time - utc_time[0] - (dates - date_i) * 86400.) / 3600.).astype(np.int)
    minutes = np.floor((utc_time - utc_time[0] - hours * 3600.0 - (dates - date_i) * 86400.) / 60.).astype(np.int)
    seconds = np.mod((utc_time - utc_time[0] - hours * 3600.0 - (dates - date_i) * 86400. - minutes * 60.0), 60.).astype(np.int)

    date_list = []
    years = np.ones(dates.__len__()) * years
    months = np.ones(dates.__len__()) * months
    for index, element in enumerate(zip(years, months, dates, hours, minutes, seconds)):
        date_list.append(str(int(element[0])) + '-' + str(int(element[1])) + '-' + str(int(element[2])) + ' ' +
                         str(int(element[3])) + ':' + str(int(element[4])) + ':' + str(element[5]))
    a_date = atime(date_list)
    xacc = dd_aird_i.xacc.compute().to_numpy()
    yacc = dd_aird_i.yacc.compute().to_numpy()
    zacc = dd_aird_i.zacc.compute().to_numpy()
    gcrs = coor.GCRS(x=xacc, y=yacc, z=zacc, representation_type='cartesian', obstime=a_date)
    itrs = gcrs.transform_to(coor.ITRS(obstime=a_date))

    itrs = np.asarray(itrs.cartesian.xyz)

    plt.style.use(['science', 'no-latex', 'high-vis'])
    fig, ax = plt.subplots(figsize=(15, 8))
    plt.plot(utc_time, itrs[0], linewidth=2, label='xacc_gcrs')
    plt.plot(utc_time, itrs[1], linewidth=2, label='yacc_gcrs')
    plt.plot(utc_time, itrs[2], linewidth=2, label='zacc_gcrs')
    plt.tick_params(labelsize=25, width=2.9)
    ax.yaxis.get_offset_text().set_fontsize(24)
    ax.xaxis.get_offset_text().set_fontsize(24)
    plt.xlabel('GPS time [s]', fontsize=20)
    plt.ylabel('$Air \quad drag \quad accelaration [km/s^2]$', fontsize=20)
    # plt.legend(fontsize=20, loc='best')
    plt.legend(fontsize=20, loc='lower left', bbox_to_anchor=(
        0, 1, 1, .1), ncol=3, mode='expand')
    plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.savefig('..//images//air_drag_itrs.png', dpi=500)
    plt.show()


if __name__ == '__main__':
    air_density(2019, 9, 16)
