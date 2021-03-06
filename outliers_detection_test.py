import datetime
import os
import time
from datetime import date
from os import read, write

import dask.dataframe as dd
import hwm93
import matplotlib.pyplot as plt
import msise00
import numba
import numpy as np
import pandas as pd

from astropy.config import configuration
from astropy import coordinates as coor
from astropy import units as u
from astropy.time import Time as atime
from astropy.coordinates import EarthLocation as el
from auromat.coordinates import transform
from matplotlib.ticker import LogLocator as LL
from matplotlib.ticker import MultipleLocator as ML
from matplotlib.ticker import ScalarFormatter as SF
from pyatmos import download_sw, nrlmsise00, read_sw
from scipy.spatial.distance import cdist


class air_drag(object):

    def __init__(self, relative_velocity, air_density, configuration, altitude):
        self.relative_velocity = relative_velocity
        self.air_density = air_density
        self.configuation = configuration
        self.altitude = altitude


@numba.jit
def diff9(f_s, time_series):
    t0 = time_series[0: -8]
    t1 = time_series[1: -7]
    t2 = time_series[2: -6]
    t3 = time_series[3: -5]
    t5 = time_series[5: -3]
    t6 = time_series[6: -2]
    t7 = time_series[7: -1]
    t8 = time_series[8:]

    diff_series = f_s / 280. * \
                  (t0 - 20. / 3. * t1 + 56. * t2 - 224. * t3 + 224. * t5 - 56. * t6 + 20. / 3. * t7 - t8)
    return diff_series


def get_run_time(func):
    def call_func(*args, **kwargs):
        begin_time = time.time()
        ret = func(*args, **kwargs)
        end_time = time.time()
        run_time = end_time - begin_time
        print(str(func.__name__) + "函数运行时间为" + str(run_time))
        return ret
    return call_func


def get_all_path(open_file_path):
    rootdir = open_file_path
    path_list = []
    list_file = os.listdir(rootdir)
    for i in range(0, len(list_file)):
        com_path = os.path.join(rootdir, list_file[i])
        # print(com_path)
        if os.path.isfile(com_path):
            if ('0860' in com_path):
                path_list.append(com_path)
    return path_list


def read_input_file_names(filepath):
    """This function is made to obtain the filenames of the input files for TAIJI-01 1D data product
        and then check the input file are continuous
    """
    # get the input files which are the coordinates of the TAIJI-01 in the earth-fixed system for a
    # whole month
    filepaths_this_month = get_all_path(filepath)

    # extract the year to check it is whether leap year or not
    year_this_month = int(filepaths_this_month[0][-25: -21])

    # extract the month to obtain the date in this month to check whether the input files are complete
    month_this_month = int(filepaths_this_month[0][-21: -19])
    month_array = check_leap_year(year_this_month, month_this_month)

    # extract the date from the input filenames
    date_this_month = np.zeros(filepaths_this_month.__len__()).astype(np.int64)
    for index_1 in np.arange(filepaths_this_month.__len__()):
        date_this_month[index_1] = int(filepaths_this_month[index_1][-19: -17])

    # check whether the date in complete or not
    if date_this_month.__len__() == month_array.__len__():
        pass
    else:
        print('Error occurred in the date of this month!')

    if month_this_month < 10:
        str_month_this_month = '0' + str(month_this_month)
    else:
        str_month_this_month = str(month_this_month)

    output_filename = '..//output//' + str_month_this_month + '//taiji-01-0860-earth-fixed-system-' + str(year_this_month) + '-' + str_month_this_month + '.txt'
    flag_filename = '..//output//' + str_month_this_month + '//taiji-01-0860-position-velocity-flag-' + str(year_this_month) + '-' + str_month_this_month + '.txt'
    air_drag_par_filename = '..//output//' + str_month_this_month + '//taiji-01-0860-air-density-' + \
        str(year_this_month) + '-' + str_month_this_month + '.txt'
    gcrs_filename = '..//output//' + str_month_this_month + '//taiji-01-0866-gcrs-' + str(year_this_month) + '-' + \
        str_month_this_month + '.txt'

    return filepaths_this_month, output_filename, flag_filename, air_drag_par_filename, gcrs_filename


def write_output_file_header(output_filename, epoch_pair, interpolated=False):
    """ This function is designed to write the header of the output file that is a imitation of the header for GRACE Follow-On date product.

    Args:
        output_filename (str): the string of the output filename.
        epoch_pair (dict[dim=n], float64): the starting epoch and the ending epoch of the series of data and other helpful messages.
    """

    local_time = time.asctime(time.localtime(time.time()))

    with open(output_filename, 'w+') as o_ef_unit:
        o_ef_unit.write('header:\n')
        o_ef_unit.write('\t' + 'dimensions: ' + str(epoch_pair['columns']) + '\n')
        o_ef_unit.write('\t' + 'global_attributes: ' + '\n')
        o_ef_unit.write('\t\t' + 'acknowledgement: ' + 'TAIJI-01 is the first grativational prospecting test satellite of China, one of whose ' +
                        'core loads is inertial sensor, therefore, the GPS data of this satellite can be used in the inversion of coefficients of earth\'s gravity field.' + '\n')
        o_ef_unit.write('\t\t' + 'creator_institution: CAS/IM\n')
        o_ef_unit.write('\t\t' + 'creator_name: TAIJI-01 data product system\n')
        o_ef_unit.write('\t\t' + 'creator_type: group\n')
        o_ef_unit.write('\t\t' + 'date_created: ' + local_time + '\n')
        o_ef_unit.write('\t\t' + 'institution: CAU\n')
        o_ef_unit.write('\t\t' + 'instrument: GPS\n')
        o_ef_unit.write('\t\t' + 'keywords: TAIJI-01, gravity field\n')
        o_ef_unit.write('\t\t' + 'processing_level: 1D\n')
        o_ef_unit.write('\t\t' + 'data_product: 0860\n')
        o_ef_unit.write('\t\t' + 'program: Taiji Programme\n')
        o_ef_unit.write('\t\t' + 'publisher_institution: CAS/ISSI\n')
        o_ef_unit.write('\t\t' + 'source: Earth-Fixed Frame trajectories for TAIJI-01\n')
        o_ef_unit.write('\t\t' + 'summary: 1-Hz trajectory states in the Earth-Fixed Frame\n')
        o_ef_unit.write('\t\t' + 'time_coverage_start: ' + str(epoch_pair['starting_epoch_gps']) + '\n')
        o_ef_unit.write('\t\t' + 'time_coverage_stop: ' + str(epoch_pair['ending_epoch_gps']) + '\n')
        o_ef_unit.write('\t\t' + 'title: TAIJI-01 Level-1D GPS Navigation Data\n')
        o_ef_unit.write('\t' + 'varaibles:\n')
        o_ef_unit.write('\t\t' + '- self_utc_time: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 1st column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Continuous seconds past ' + str(epoch_pair['starting_epoch_gps']) + '\n')
        o_ef_unit.write('\t\t\t' + 'unit: second\n')
        o_ef_unit.write('\t\t' + '- self_rcv_time: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 2nd column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Continuous seconds past ' + epoch_pair['starting_epoch_rcv'] + '\n')
        o_ef_unit.write('\t\t\t' + 'unit: second\n')
        o_ef_unit.write('\t\t' + '- xpos: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 3rd column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Position, X value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km\n')
        o_ef_unit.write('\t\t' + '- ypos: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 4th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Position, Y value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km\n')
        o_ef_unit.write('\t\t' + '- zpos: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 5th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Position, Z value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km\n')
        o_ef_unit.write('\t\t' + '- xvel: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 6th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Velocity, X value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s\n')
        o_ef_unit.write('\t\t' + '- yvel: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 7th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Velocity, Y value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s\n')
        o_ef_unit.write('\t\t' + '- zvel: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 8th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Velocity, Z value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s\n')
        o_ef_unit.write('\t\t' + '- xacc: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 9th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Accelaration, X value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s^2\n')
        o_ef_unit.write('\t\t' + '- yacc: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 10th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Accelaration, Y value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s^2\n')
        o_ef_unit.write('\t\t' + '- zacc: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 11th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Accelaration, Z value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s^2\n')
        o_ef_unit.write('\t\t' + '- vac_qualflg: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 12th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: the symbol which indicates that missing several epoches before this epoch\n')
        o_ef_unit.write('\t\t\t' + 'interpretation: 0 is the rightmost\n')
        o_ef_unit.write('\t\t' + '- xvel_diffed: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 13th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Velocity diffed from position x, X value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s\n')
        o_ef_unit.write('\t\t' + '- yvel_diffed: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 14th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Velocity diffed from position y, Y value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s\n')
        o_ef_unit.write('\t\t' + '- zvel_diffed: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 15th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: Velocity diffed from position z, Z value\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s\n')
        o_ef_unit.write('\t\t' + '- pj_qualflg: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 16th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: the symbol which indicates position jump\n')
        o_ef_unit.write('\t\t\t' + 'interpretation: 0 is the rightmost\n')
        o_ef_unit.write('\t\t' + '- vj_qualflg: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 17th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: the symbol which indicates velocity jump\n')
        o_ef_unit.write('\t\t\t' + 'interpretation: 0 is the rightmost\n')
        o_ef_unit.write('# End of header\n')


def creat_qualflg(df, flag_filename):
    """ This function is designed to creat qualflg for every time epoch

    Args:
        df (dataframe): this is the dataframe containing the time, positions and velocities.
    """

    utc_time = df['utc_time'].compute().to_numpy()
    utc_time = utc_time - utc_time[0]
    position_incr = np.zeros([df.__len__(), 3]); velocity_incr = np.zeros([df.__len__(), 3])
    trajectory = df[['xpos', 'ypos', 'zpos', 'xvel', 'yvel', 'zvel']].compute().to_numpy()
    trajectory_flag = np.abs(np.diff(np.diff(trajectory, axis=0), axis=0))
    trajectory_flag_unabs = np.diff(np.diff(trajectory, axis=0), axis=0)
    position_flag = np.logical_or(trajectory_flag[:, 0] > 1, trajectory_flag[:, 1] > 1, trajectory_flag[:, 2] > 1).astype(np.int64)
    velocity_flag = np.logical_or(trajectory_flag[:, 3] > 0.0025, trajectory_flag[:, 4] > 0.0025, trajectory_flag[:, 5] > 0.0025).astype(np.int64)
    position_mask = np.where(position_flag == 1)[0]
    position_mask = position_mask[np.arange(1, position_mask.__len__(), 2)]
    position_flag[position_mask] = 0
    position_incr[position_mask + 1, 0] = trajectory_flag_unabs[position_mask - 1, 0]
    position_incr[position_mask + 1, 1] = trajectory_flag_unabs[position_mask - 1, 1]
    position_incr[position_mask + 1, 2] = trajectory_flag_unabs[position_mask - 1, 2]
    velocity_mask = np.where(velocity_flag == 1)[0]
    velocity_mask = velocity_mask[np.arange(1, velocity_mask.__len__(), 2)]
    velocity_flag[velocity_mask] = 0
    velocity_incr[velocity_mask + 1, 0] = trajectory_flag_unabs[velocity_mask - 1, 3]
    velocity_incr[velocity_mask + 1, 1] = trajectory_flag_unabs[velocity_mask - 1, 4]
    velocity_incr[velocity_mask + 1, 2] = trajectory_flag_unabs[velocity_mask - 1, 5]
    position_flag = np.r_[0, 0, position_flag]; velocity_flag = np.r_[0, 0, velocity_flag]

    # append new column
    if position_flag.__len__() != df.__len__():
        raise ValueError('Error appending the position flag to the dataframe.')
    else:
        # dd_gps['position_flag'] = dd.from_array(position_flag)
        # dd_gps['velocity_flag'] = dd.from_array(velocity_flag)
        np.savetxt(flag_filename, np.c_[position_flag, velocity_flag, position_incr, velocity_incr], delimiter=' ', newline='\n',
                   header='position velocity posotion_increment velocity_increment', fmt='%1d %1d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f')

    # plt.style.use(['science', 'no-latex', 'high-vis'])
    # fig, ax1 = plt.subplots(figsize=(15, 8))
    # ax1.plot(trajectory_flag[:, 0], linewidth=2, label='xpos_dot', marker='o')
    # ax1.plot(trajectory_flag[:, 1], linewidth=2, label='ypos_dot')
    # ax1.plot(trajectory_flag[:, 2], linewidth=2, label='zpos_dot')
    # ax1.tick_params(axis='x', which='minor', rotation=45, labelsize=15)
    # ax1.yaxis.get_offset_text().set_fontsize(24)
    # ax1.xaxis.get_offset_text().set_fontsize(24)
    # ax1.set_xlabel('time [s]', fontsize=24)
    # ax1.set_ylabel('velocity [km/s]', fontsize=24)
    # # plt.legend(fontsize=20, loc='lower left', bbox_to_anchor=(0, 1, 1, .1), ncol=4, mode='expand')
    # plt.legend(fontsize=20, loc='best')
    # plt.tick_params(labelsize=25, axis='both')
    # plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    # plt.gca().spines['left'].set_linewidth(2)
    # plt.gca().spines['top'].set_linewidth(2)
    # plt.gca().spines['right'].set_linewidth(2)
    # plt.gca().spines['bottom'].set_linewidth(2)
    #
    # fig, ax1 = plt.subplots(figsize=(15, 8))
    # ax1.plot(trajectory_flag[:, 3], linewidth=2, label='xvel_dot', marker='o')
    # ax1.plot(trajectory_flag[:, 4], linewidth=2, label='yvel_dot')
    # ax1.plot(trajectory_flag[:, 5], linewidth=2, label='zvel_dot')
    # ax1.tick_params(axis='x', which='minor', rotation=45, labelsize=15)
    # ax1.yaxis.get_offset_text().set_fontsize(24)
    # ax1.xaxis.get_offset_text().set_fontsize(24)
    # ax1.set_xlabel('time [s]', fontsize=24)
    # ax1.set_ylabel('velocity [km/s]', fontsize=24)
    # # plt.legend(fontsize=20, loc='lower left', bbox_to_anchor=(0, 1, 1, .1), ncol=4, mode='expand')
    # plt.legend(fontsize=20, loc='best')
    # plt.tick_params(labelsize=25, axis='both')
    # plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    # plt.gca().spines['left'].set_linewidth(2)
    # plt.gca().spines['top'].set_linewidth(2)
    # plt.gca().spines['right'].set_linewidth(2)
    # plt.gca().spines['bottom'].set_linewidth(2)

    # plt.show()
    # return position_flag, velocity_flag


def air_drag_par(df, filename, sw_data, gcrs_filename, date_i):
    """ This function is designed to calculate the wind speed for the satellite.

    Args:
        df (dask.datafram): dataframe containing different information
        filename (str): the output file name
    """
    # extract year month date hour minute second
    year = np.int(filename[-11: -7])
    year = np.dot(np.ones(df.__len__()), year).astype(np.int)
    month = np.int(filename[-6: -4])
    month = np.dot(np.ones(df.__len__()), month).astype(np.int)
    utc_time = df.utc_time.compute().to_numpy()
    xpos = df.xpos.compute().to_numpy()
    ypos = df.ypos.compute().to_numpy()
    zpos = df.zpos.compute().to_numpy()
    date = np.floor((utc_time - utc_time[0]) / 86400.).astype(np.int) + date_i
    hour = np.floor(
        (utc_time - utc_time[0] - (date - date_i) * 86400.) / 3600.).astype(np.int)
    minute = np.floor(
        (utc_time - utc_time[0] - hour * 3600.0 - (date - date_i) * 86400.) / 60.).astype(np.int)
    second = np.mod((utc_time - utc_time[0] - hour * 3600.0 - (
        date - date_i) * 86400. - minute * 60.0), 60.).astype(np.int)
    # using astropy transform ecef to geodetic
    pos_el = el.from_geocentric(x=xpos, y=ypos, z=zpos, unit='km')
    lla_el = pos_el.to_geodetic()
    lat = np.asarray(lla_el.lat)
    lon = np.asarray(lla_el.lon)
    alt = np.asarray(lla_el.height)

    # wind speed
    winds = np.zeros([df.__len__(), 2])
    air_density = np.zeros(df.__len__())
    for i, y in enumerate(year):
        y = np.int(y); h = np.int(hour[i])
        mon = np.int(month[i]); d = np.int(date[i]); minu = np.int(minute[i]); sec = np.int(second[i])
        date_time = datetime.datetime(y, mon, d, h, minu, sec)
        para_input, para_output = nrlmsise00(date_time, lat[i], lon[i], alt[i], sw_data)
        air_density[i] = para_output['Density']['RHO[kg/m^3]']  # unit = km / m^(-3)
        reg= hwm93.run(
            date_time, altkm=alt[i], glat=lat[i], glon=lon[i], f107=para_input['f107Daily[10^-22 W/m^2/Hz]'],
            f107a=para_input['f107Average[10^-22 W/m^2/Hz]'], ap=para_input['ApDaily'])
        winds[i, 0] = reg.meridional.values
        winds[i, 1] = reg.zonal.values


    # cartesian 2 spherical
    slat, slon = transform.cartesian_to_spherical(xpos, ypos, zpos, with_radius=False)
    # lat
    ns_morethan_0 = np.where(winds[:, 0] >= 0)
    ns_lessthan_0 = np.where(winds[:, 0] < 0)
    slat[ns_morethan_0] = np.abs(slat[ns_morethan_0] - np.pi / 2.0)
    slat[ns_lessthan_0] = -np.abs(slat[ns_lessthan_0] - np.pi / 2.0)
    ns_x, ns_y, ns_z = transform.spherical_to_cartesian(winds[:, 0], slat, slon)
    # lon
    lon_lessthan_0 = np.where(slon < 0)
    slon[lon_lessthan_0] = 2 * np.pi + slon[lon_lessthan_0]
    slon[ns_morethan_0] = slon[ns_morethan_0] + np.pi / 2.0
    slon[ns_lessthan_0] = slon[ns_lessthan_0] - np.pi / 2.0
    lon_lessthan_0 = np.where(slon < 0)
    slon[lon_lessthan_0] = 2 * np.pi + slon[lon_lessthan_0]
    lon_morethan_2 = np.where(slon >= 2 * np.pi)
    slon[lon_morethan_2] = slon[lon_morethan_2] - 2 * np.pi
    lon_morethan_p = np.where(slon > np.pi)
    slon[lon_morethan_p] = slon[lon_morethan_p] - 2 * np.pi
    ew_x, ew_y, ew_z = transform.spherical_to_cartesian(winds[:, 1], np.zeros(slon.__len__()), slon)
    winds_x = ew_x + ns_x; winds_y = ew_y + ns_y; winds_z = ew_z + ns_z
    winds_vector = np.c_[winds_x, winds_y, winds_z]
    angular_vector = np.asarray([0.0, 0.0, 7.2921158553e-5])

    spos = np.c_[xpos, ypos, zpos]
    svel = np.c_[df.xvel.compute().to_numpy(), df.yvel.compute().to_numpy(), df.zvel.compute().to_numpy()]
    relative_velocity = (svel * 1000. - winds_vector - np.cross(angular_vector, spos * 1000.)) / 1000.
    itrs2gcrs(year, month, date, hour, minute, second, spos, relative_velocity, utc_time, gcrs_filename)

    return relative_velocity, air_density


def init_pyatmos():
    """ This is a function to initiate the package-pyatmos

    Returns:
        [type]: [description]
    """
    swfile = download_sw()
    swdata = read_sw(swfile)

    return swdata


def itrs2gcrs(years, months, dates, hours, minutes, seconds, pos, vel, utc_time, gcrs_filename):
    """ This function is to transform the coordinates and velocities from itrs into gcrs
    """
    date_list = []
    for index, element in enumerate(zip(years, months, dates, hours, minutes, seconds)):
        date_list.append(str(int(element[0])) + '-' + str(int(element[1])) + '-' + str(int(element[2])) + ' ' + \
            str(int(element[3])) + ':' + str(int(element[4])) + ':' + str(element[5]))

    a_date = atime(date_list)
    pos = pos * u.km; vel = vel * u.km
    itrs1 = coor.ITRS(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2], representation_type='cartesian', obstime=a_date)
    gcrs1 = itrs1.transform_to(coor.GCRS(obstime=a_date))
    itrs2 = coor.ITRS(x=vel[:, 0], y=vel[:, 1], z=vel[:, 2], representation_type='cartesian', obstime=a_date)
    gcrs2 = itrs2.transform_to(coor.GCRS(obstime=a_date))

    gcrs_xyz1 = np.asarray(gcrs1.cartesian.xyz)
    gcrs_xyz2 = np.asarray(gcrs2.cartesian.xyz)

    dd_igrf = dd.from_pandas(pd.DataFrame({'utc_time': utc_time, 'gcrs_x': gcrs_xyz1[0],
                                           'gcrs_y': gcrs_xyz1[1], 'gcrs_z': gcrs_xyz1[2],
                                           'gcrs_rvx': gcrs_xyz2[0],
                                           'gcrs_rvy': gcrs_xyz2[1],
                                           'gcrs_rvz': gcrs_xyz2[2]}),
                             npartitions=1)
    write_gcrs_file_header(gcrs_filename, dd_igrf)
    dd_igrf.to_csv(gcrs_filename, single_file=True, sep='\t', index=False, mode='a+')


def write_gcrs_file_header(output_filename, df):
    """ This function is designed to write the header of the output file that is a imitation of the header for GRACE Follow-On date product.

    Args:
        output_filename (str): the string of the output filename.
        epoch_pair (dict[dim=n], float64): the starting epoch and the ending epoch of the series of data and other helpful messages.
    """

    local_time = time.asctime(time.localtime(time.time()))

    with open(output_filename, 'w+') as o_ef_unit:
        o_ef_unit.write('header:\n')
        o_ef_unit.write('\t' + 'dimensions: ' +
                        str(df.__len__()) + '\n')
        o_ef_unit.write('\t' + 'global_attributes: ' + '\n')
        o_ef_unit.write('\t\t' + 'acknowledgement: ' + 'TAIJI-01 is the first grativational prospecting test satellite of China, one of whose ' +
                        'core loads is inertial sensor, therefore, the Beidou data of this satellite can be used in the inversion of coefficients of earth\'s gravity field.' + '\n')
        o_ef_unit.write('\t\t' + 'creator_institution: CAS/IM\n')
        o_ef_unit.write(
            '\t\t' + 'creator_name: TAIJI-01 data product system\n')
        o_ef_unit.write('\t\t' + 'creator_type: group\n')
        o_ef_unit.write('\t\t' + 'date_created: ' + local_time + '\n')
        o_ef_unit.write('\t\t' + 'institution: CAU\n')
        o_ef_unit.write('\t\t' + 'instrument: Beidou receiver\n')
        o_ef_unit.write('\t\t' + 'keywords: TAIJI-01, gravity field\n')
        o_ef_unit.write('\t\t' + 'processing_level: 1D\n')
        o_ef_unit.write('\t\t' + 'data_product: 0860\n')
        o_ef_unit.write('\t\t' + 'program: Taiji Programme\n')
        o_ef_unit.write('\t\t' + 'publisher_institution: CAS/ISSI\n')
        o_ef_unit.write(
            '\t\t' + 'source: Earth-Fixed Frame trajectories for TAIJI-01\n')
        o_ef_unit.write(
            '\t\t' + 'summary: 1-Hz trajectory states in the Earth-Fixed Frame\n')
        o_ef_unit.write('\t\t' + 'time_coverage_start: ' +
                        str(df.utc_time.compute().to_numpy()[0]) + '\n')
        o_ef_unit.write('\t\t' + 'time_coverage_stop: ' +
                        str(df.utc_time.compute().to_numpy()[-1]) + '\n')
        o_ef_unit.write('\t\t' + 'title: TAIJI-01 Level-1D Beidou Navigation Data\n')
        o_ef_unit.write('\t' + 'varaibles:\n')
        o_ef_unit.write('\t\t' + '- self_utc_time: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 1st column\n')
        o_ef_unit.write('\t\t\t' + 'unit: second\n')
        o_ef_unit.write('\t\t' + '- xpos: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 2nd column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: the x position in gcrs\n')
        o_ef_unit.write('\t\t\t' + 'unit: km\n')
        o_ef_unit.write('\t\t' + '- ypos: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 3rd column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: the y position in gcrs\n')
        o_ef_unit.write('\t\t\t' + 'unit: km\n')
        o_ef_unit.write('\t\t' + '- zpos: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 4th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: the z position in gcrs\n')
        o_ef_unit.write('\t\t\t' + 'unit: km\n')
        o_ef_unit.write('\t\t' + '- xrvel: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 5th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: the relative velocity x in gcrs\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s\n')
        o_ef_unit.write('\t\t' + '- yrvel: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 6th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: the relative velocity y in gcrs\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s\n')
        o_ef_unit.write('\t\t' + '- zrvel: \n')
        o_ef_unit.write('\t\t\t' + 'comment: 7th column\n')
        o_ef_unit.write('\t\t\t' + 'long_name: the relative velocity z in gcrs\n')
        o_ef_unit.write('\t\t\t' + 'unit: km/s\n')
        o_ef_unit.write('# End of header\n')


@get_run_time
def outliers_detection():
    """ This function is designed to find the outliers for the GPS data in earth-fixed system of Taiji-01.
    """
    date_init = 16.0

    # get the input file names and the output file name
    filepaths_this_month, output_filename, flag_filename, air_drag_par_filename, gcrs_filename = read_input_file_names('..//input//09')

    # read the input files all together
    dd_gps = dd.read_csv(urlpath=filepaths_this_month, sep=',', header=None,
                         engine='c', skiprows=2, storage_options=dict(auto_mkdir=False), names=['rcv_time', 'gn001',
                         'gn002', 'gn003', 'gn004', 'gn005', 'gn006', 'gn007', 'gn008', 'gn009', 'gn010', 'gn011',
                         'gn012', 'gn013', 'gn014', 'gn015-0', 'gn15-1', 'gn16', 'gn17', 'gn18', 'gn19', 'gn20', 'gn21',
                         'gn022', 'gn023', 'gn024', 'gn025', 'gn026', 'gn027', 'gn028', 'gn029', 'gn030', 'gn031', 'gn032',
                         'gn033', 'gn034', 'utc_time', 'xpos', 'ypos', 'zpos', 'xvel', 'yvel', 'zvel', 'gn042', 'gn043',
                         'gn044', 'gn045', 'gn046', 'gn047', 'gn048', 'gn049', 'gn050', 'gn051', 'gn052', 'gn053', 'gn054',
                         'gn055', 'gn056', 'gn057', 'WTF1', 'WTF2', 'WTF3'], dtype={'xpos': np.float64, 'ypos': np.float64, 'zpos': np.float64,
                         'xvel': np.float64, 'yvel': np.float64, 'zvel': np.float64, 'gn056': 'object'}, encoding='gb2312')

    # compute the self time frame
    dd_gps = dd_gps.drop_duplicates(subset=['rcv_time'])
    dd_gps = dd_gps.drop_duplicates(subset=['utc_time'])
    dd_gps = dd_gps[['rcv_time', 'xpos', 'ypos', 'zpos', 'xvel', 'yvel', 'zvel', 'utc_time']]
    dd_gps['T-bfore_rcv'] = dd_gps['rcv_time'].str.partition('T')[0]
    dd_gps['T-after_rcv'] = dd_gps['rcv_time'].str.partition('T')[2]
    dd_gps['year_rcv'] = dd_gps['T-bfore_rcv'].str.partition('-')[0].astype(np.int32)
    dd_gps['month_rcv'] = dd_gps['T-bfore_rcv'].str.partition('-')[2].str.partition('-')[0].astype(np.int32)
    dd_gps['date_rcv'] = dd_gps['T-bfore_rcv'].str.partition('-')[2].str.partition('-')[2].astype(np.int32)
    dd_gps['hour_rcv'] = dd_gps['T-after_rcv'].str.partition(':')[0].astype(np.int32)
    dd_gps['minute_rcv'] = dd_gps['T-after_rcv'].str.partition(':')[2].str.partition(':')[0].astype(np.int32)
    dd_gps['second_rcv'] = dd_gps['T-after_rcv'].str.partition(':')[2].str.partition(':')[2].astype(np.float64)
    dd_gps['self_rcv_time'] = (dd_gps['date_rcv'] - 1.0) * 86400.0 + dd_gps['hour_rcv'] * 3600.0 + dd_gps['minute_rcv'] * 60.0 + dd_gps['second_rcv']

    # sort the dataframe with the gps time
    dd_gps = dd_gps[dd_gps['utc_time'] > 300000.0]
    dd_gps = dd_gps.nsmallest(dd_gps.__len__(), 'utc_time')

    epoch_pair = {'starting_epoch_gps': dd_gps['utc_time'].compute().to_numpy()[0], 'ending_epoch_gps': dd_gps['utc_time'].compute().to_numpy()[-1],
                  'columns': dd_gps.__len__(), 'instrument': 'Beidou', 'data_type': 'earth_fixed_system_positions_and_velocities',
                  'version': 1, 'starting_epoch_rcv': dd_gps['rcv_time'].compute().to_numpy()[0], 'ending_epoch_rcv': dd_gps['rcv_time'].compute().to_numpy()[-1]}

    # calculate the wind speed
    sw_data = init_pyatmos()
    relative_velocity, air_density = air_drag_par(dd_gps, air_drag_par_filename, sw_data, gcrs_filename, date_init)
    np.savetxt(air_drag_par_filename, air_density, fmt='%.10e', delimiter=' ', newline='\n',
               header='Air density (kg / m^3)')

    # creat the qualflg
    creat_qualflg(dd_gps, flag_filename)

    # write the output file
    write_output_file_header(output_filename=output_filename, epoch_pair=epoch_pair)

    dd_gps.to_csv(output_filename, single_file=True, sep='\t', index=False, mode='a+', columns=[
        'utc_time', 'self_rcv_time', 'xpos', 'ypos', 'zpos', 'xvel', 'yvel', 'zvel'])


def check_leap_year(year, month):
    """This function is designed to check the input year is wether leap year or not. If the input
        year is leap year, the values of February will be add with a number------29;

    Args:
        year (int64): the year of the input filename.
        month (int64): the month of the input filename.
    """
    month_dict = {
        1: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
        2: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28],
        3: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
        4: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30],
        5: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
        6: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30],
        7: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
        8: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
        9: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30],
        10: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31],
        11: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30],
        12: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
    }

    if (year % 4) == 0:
        if (year % 100) == 0:
            if (year % 400) == 0:
                month_dict[2].append(29)
        else:
            month_dict[2].append(29)
    return month_dict[month]


if __name__ == '__main__':
    outliers_detection()
