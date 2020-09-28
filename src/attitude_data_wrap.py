import os
import time
import numba

import dask.dataframe as dd
import numpy as np
import pandas as pd
import spiceypy as sp
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.signal import kaiserord, firwin, filtfilt, lfilter, freqz, welch



class Att(object):
    # define a class named att
    def __init__(self, df, days):
        self.df = df
        self.days = days


    def date_reform(self, form):
        """ This is a function to reform the date from string form to yyyy|mm|dd|HH|MM|SS

        Args:
            form (string): the form of the input date string

        Returns:
            float64: yyyy|mm|dd|HH|MM|SS
        """
        if form != 'yyyy-mm-ddTHH:MM:SS':
            raise ValueError('The input time form can not be recognised.')
        else:
            self.df['T-bfore_rcv'] = self.df['rcv_time'].str.partition('T')[0]
            self.df['T-after_rcv'] = self.df['rcv_time'].str.partition('T')[2]
            self.df['year_rcv'] = self.df['T-bfore_rcv'].str.partition('-')[
                0].astype(np.int32)
            self.df['month_rcv'] = self.df['T-bfore_rcv'].str.partition(
                '-')[2].str.partition('-')[0].astype(np.int32)
            self.df['date_rcv'] = self.df['T-bfore_rcv'].str.partition(
                '-')[2].str.partition('-')[2].astype(np.int32)
            self.df['hour_rcv'] = self.df['T-after_rcv'].str.partition(':')[
                0].astype(np.int32)
            self.df['minute_rcv'] = self.df['T-after_rcv'].str.partition(
                ':')[2].str.partition(':')[0].astype(np.int32)
            self.df['second_rcv'] = self.df['T-after_rcv'].str.partition(
                ':')[2].str.partition(':')[2].astype(np.float64)
            self.df['self_rcv_time'] = (self.df['date_rcv'] - 1.0) * 86400.0 + self.df['hour_rcv'] * \
                3600.0 + self.df['minute_rcv'] * 60.0 + self.df['second_rcv']

            # sort the dataframe with the gps time
            # self.df = self.df[self.df['gps_time'] > 300000.0]
            self.df['index'] = self.df['self_rcv_time']
            self.df = self.df.set_index('index')
            self.df = self.df.nsmallest(self.df.__len__(), 'self_rcv_time')


    def equidistant_quantity(self, flags, interval):
        """ This function is designed to interpolate the original gps_time physical quantity sereis
            to an equidistant one

        Args:
            flags: the column(s) which are interpolating
            interval: the interval between two near epoches

        Raises:
            ValueError: [description]

        Return:
            df: new dataframe containing the columns after interpolation. Note that the length of
                this new dataframe is not necessarily equal to the length of the input one
        """
        # generate the equidistant series
        floor_init_epoch = np.floor(self.df.self_rcv_time.compute().to_numpy()[0]) + 52588800.0
        time_span = np.arange(floor_init_epoch, floor_init_epoch + 86400. * self.days,
                              interval).astype(np.float64)

        flags_list = []
        temp = np.zeros([time_span.__len__(), 7])
        temp[:, 0] = time_span

        for index, flag in enumerate(flags):
            interp_func = interpolate.interp1d(self.df['gps_time'].compute().to_numpy(),
                                               self.df[flag].compute().to_numpy(), kind='linear',
                                               bounds_error=False, fill_value='extrapolate')
            flags_list.append(flag + '_interp')
            temp[:, index + 1] = interp_func(time_span)

        self.df = dd.from_pandas(pd.DataFrame({'gps_time': temp[:, 0],
                                               flags_list[0]: temp[:, 1],
                                               flags_list[1]: temp[:, 2],
                                               flags_list[2]: temp[:, 3],
                                               flags_list[3]: temp[:, 4],
                                               flags_list[4]: temp[:, 5],
                                               flags_list[5]: temp[:, 6]}), npartitions=1)

        return flags_list


    def down_sample(self, flags, cutoff_hz, fq, ripple_db=300.):
        """ This method is designed to down sample some time series

        Args:
            fq (float64): niquist frequency
        """
        # the desired width of the transition from pass to stop
        width = 0.52 / fq

        # the desired attenuation in the stop band in db: ripple_db

        # compute the kaiser parameter for the fir filter
        n, beta = kaiserord(ripple_db, width)
        print('The length of the lowpass filter is ', n, '.')

        # use firwin with a kaiser window
        taps = firwin(n, cutoff_hz, window=('kaiser', beta),
                      pass_zero='lowpass', nyq=fq)

        # use filtfilt to filter x with the fir filter
        filtered_list = []
        for index, flag in enumerate(flags):
            filtered_x = filtfilt(taps, 1.0, self.df[flag].compute().to_numpy())
            print(filtered_x.__len__())
            self.df = self.df.reset_index().set_index('index')
            filtered_list.append(flag + '_filt')
            self.df[filtered_list[index]] = 0
            self.df[filtered_list[index]] = self.df[filtered_list[index]].compute() + filtered_x

        # plot the frequency response
        plot_freq_response(taps, fq)
        self.df = self.df.nsmallest(self.df.__len__(), 'gps_time')

        return filtered_list


    def write_file(self, flags, outputfilename):
        """ This function is designed to write the desired columns into files

        Args:
            flags (string list): desired columns

        Raises:
            ValueError: [description]

        Returns:
            [type]: [description]
        """
        flags.append('gps_time')
        flags[0], flags[1: ] = flags[-1], flags[0: -1]

        print(self.df[flags[1]].compute().to_numpy()[-1])
        self.df.to_csv(outputfilename, single_file=True, sep='\t', index=False, mode='a+',
                       columns=flags)


    def write_output_file_header(self, output_filename):
        """ This function is designed to write the header of the output file that is a imitation of the header for GRACE Follow-On date product.

        Args:
            output_filename (str): the string of the output filename.
            epoch_pair (dict[dim=n], float64): the starting epoch and the ending epoch of the series of data and other helpful messages.
        """

        local_time = time.asctime(time.localtime(time.time()))

        with open(output_filename, 'w+') as o_ef_unit:
            o_ef_unit.write('header:\n')
            o_ef_unit.write('\t' + 'dimensions: ' + str(self.df.__len__()) + '\n')
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
            o_ef_unit.write('\t\t' + 'time_coverage_start: ' + str(self.df.gps_time.compute().to_numpy()[0]) + '\n')
            o_ef_unit.write('\t\t' + 'time_coverage_stop: ' +
                            str(self.df.gps_time.compute().to_numpy()[0]) + '\n')
            o_ef_unit.write('\t\t' + 'title: TAIJI-01 Level-1D GPS Navigation Data\n')
            o_ef_unit.write('\t' + 'varaibles:\n')
            o_ef_unit.write('\t\t' + '- self_gps_time: \n')
            o_ef_unit.write('\t\t\t' + 'comment: 1st column\n')
            o_ef_unit.write('\t\t\t' + 'unit: second\n')
            o_ef_unit.write('\t\t' + '- igrf_eul_x: \n')
            o_ef_unit.write('\t\t\t' + 'comment: 2nd column\n')
            o_ef_unit.write('\t\t\t' + 'unit: degree\n')
            o_ef_unit.write('\t\t' + '- igrf_eul_y: \n')
            o_ef_unit.write('\t\t\t' + 'comment: 3rd column\n')
            o_ef_unit.write('\t\t\t' + 'unit: degree\n')
            o_ef_unit.write('\t\t' + '- igrf_eul_z: \n')
            o_ef_unit.write('\t\t\t' + 'comment: 4th column\n')
            o_ef_unit.write('\t\t\t' + 'unit: degree\n')
            o_ef_unit.write('\t\t' + '- tf_eul_x: \n')
            o_ef_unit.write('\t\t\t' + 'comment: 5th column\n')
            o_ef_unit.write('\t\t\t' + 'unit: degree\n')
            o_ef_unit.write('\t\t' + '- tf_eul_y: \n')
            o_ef_unit.write('\t\t\t' + 'comment: 6th column\n')
            o_ef_unit.write('\t\t\t' + 'unit: degree\n')
            o_ef_unit.write('\t\t' + '- tf_eul_z: \n')
            o_ef_unit.write('\t\t\t' + 'comment: 7th column\n')
            o_ef_unit.write('\t\t\t' + 'unit: degree\n')
            o_ef_unit.write('# End of header\n')


class InputFile(object):
    # define a class named InputFile
    def __init__(self, flag, year, month, dir_path, df, path_list, output_filename, days):
        self.flag = flag
        self.month = month
        self.dir_path = dir_path
        self.df = df
        self.path_list = path_list
        self.year = year
        self.output_filename = output_filename
        self.days = days


    # return all the attitude file paths
    def io_file_paths(self):
        """ This is a function to return the file paths for the flag data.

        Args:
            flag (str): the flag indicating file type

        Return:
            path_list (string array):
        """
        if self.month < 10:
            month_str = '0' + str(self.month)
        else:
            month_str = str(self.month)

        rootdir = self.dir_path + month_str

        self.path_list = []
        list_file = os.listdir(rootdir)
        for i in range(0, len(list_file)):
            com_path = os.path.join(rootdir, list_file[i])
            if os.path.isfile(com_path):
                if (self.flag in com_path):
                    if 'swp' not in com_path: # prevent vim openness intercepting
                        self.path_list.append(com_path)

        self.days = self.path_list.__len__()

        try:
            # extract the year of this directory to check it's leap year or not
            year_this_month = int(self.path_list[0][-25: -21])
            # extract the month to obtain the date in this month to check whether the input files are complete
            month_this_month = int(self.path_list[0][-21: -19])
            month_array = check_leap_year(year_this_month, month_this_month)

            # extract the date from the input filenames
            date_this_month = np.zeros(self.path_list.__len__()).astype(np.int64)
            for index_1 in np.arange(self.path_list.__len__()):
                date_this_month[index_1] = int(self.path_list[index_1][-19: -17])

            # check whether the date in complete or not
            if date_this_month.__len__() != month_array.__len__():
                print('Error occurred in the date of this month!')

            # output file name
            self.output_filename = '..//output//taiji-01-0811-attitude-' + \
                str(year_this_month) + '-' + month_str + '.txt'

            self.df = self.df.nsmallest(self.df.__len__(), 'self_rcv_time')
        except:
            raise ValueError('The format of the input file(s) is not correct')

        return self.path_list, self.output_filename, self.days


    # read attitude data from files
    def read_data(self, flag):
        """ This is a function to read attitude data from the files. Due to the vast differences for
            different file types, the read_file_type_data function is varied from type to type

        Args:
            filepaths (string array): the string array containing the attitude file paths
        """

        # hedaer dictionary for different files
        header_dict = {'0811': ['卫星时间', '无拖曳与姿态控制工作模式', '姿态确定模式', '角速度确定模式',
                       '轨道有效标志', '保留1', '陀螺已连续积分时间(s)', '保留2', '冷推控制模式', '量轮控制模式',
                       '磁力矩器控制模式', '自主转入对地定向使能状态', '冷气微推力器速率阻尼允许/禁止状态',
                       '星敏仲裁屏蔽状态', 'Wbi 磁控速率阻尼允许/禁止状态', '双星敏定姿允许/禁止状态',
                       '保留3', '保留4', '星敏感器角速度工作状态字', '进安全模式原因标识',
                       '用于控制的本体系太阳矢量状态','用于控制的本体系地磁矢量状态', '星敏感器四元数工作状态字',
                       '光纤陀螺工作状态字', '太敏工作状态字', '磁强计工作状态字', '保留5', '保留6', '保留7', '保留8',
                       '磁力矩器输出屏蔽状态', '磁力矩器 X 轴极性', '磁力矩器 Y 轴极性', '磁力矩器 Z 轴极性',
                       '磁力矩器发出指令电压 X大小', '磁力矩器发出指令电压 Y大小', '磁力矩器发出指令电压 Z大小',
                       '磁力矩器发出指令电压 X方向', '磁力矩器发出指令电压 Y方向', '磁力矩器发出指令电压 Z方向',
                       '动量轮工作状态', '动量轮转速有效标志', '动量轮通讯状态标志', '星敏 A 接口状态字',
                       '星敏 A 单机状态字', '星敏 B 接口状态字', '星敏 B 单机状态字', '陀螺 A 接口状态字',
                       '陀螺 A 单机状态字', '陀螺 B 接口状态字', '陀螺 B 单机状态字', '冷气微推力器 1 工作状态',
                       '冷气微推力器 2 工作状态', '冷气微推力器 3 工作状态', '冷气微推力器 4 工作状态',
                       '冷气微推力器 5 工作状态', '冷气微推力器 6 工作状态', '冷气微推力器 7 工作状态',
                       '冷气微推力器 8 工作状态', '冷气微推力器输出屏蔽状态', '冷气微推力器 1 指令喷气时长(ms)',
                       '冷气微推力器 2 指令喷气时长(ms)', '冷气微推力器 3 指令喷气时长(ms)',
                       '冷气微推力器 4 指令喷气时长(ms)', '冷气微推力器 5 指令喷气时长(ms)',
                       '冷气微推力器 6 指令喷气时长(ms)', '冷气微推力器 7 指令喷气时长(ms)',
                       '冷气微推力器 8 指令喷气时长(ms)', '太阳矢量与-Yb 轴的夹角', '轨道系姿态角 X',
                       '轨道系姿态角 Y', '轨道系姿态角 Z', '惯性系姿态角速度 X', '惯性系姿态角速度 Y',
                       '惯性系姿态角速度 Z', '轨道系角速度 Wbox', '轨道系角速度 Wboy', '轨道系角速度 Wboz',
                       '导引律控制坐标系欧拉角X', '导引律控制坐标系欧拉角Y', '导引律控制坐标系欧拉角Z',
                       '惯性系姿态角 X', '惯性系姿态角 Y', '惯性系姿态角 Z', '本体系磁矢量 X', '本体系磁矢量 Y',
                       '本体系磁矢量 Z','本体系太阳矢量 X','本体系太阳矢量 Y', '本体系太阳矢量 Z', 'WTF',
                       'WTF1', 'WTF2'],
                       '0111': ['rcv_time', 'eps_time']}

        # the dtype for different columns
        dtype_dict = {'0811': {'卫星时间': str, '惯性系姿态角 X': np.float64,
                               '惯性系姿态角 Y': np.float64, '惯性系姿态角 Z': np.float64,
                               '轨道系姿态角 X': np.float64, '轨道系姿态角 Y': np.float64,
                               '轨道系姿态角 Z': np.float64},
                      '0111': {'rcv_time': np.float64, 'eps_time': np.float64}}

        # the new column
        new_header_dict = {'0811': ['rcv_time', 'igrf_eul_x', 'igrf_eul_y', 'igrf_eul_z',
                                    'tf_eul_x', 'tf_eul_y', 'tf_eul_z'],
                           '0111': ['rcv_time', 'eps_time']}

        # used columns
        used_column_dict = {'0811': ['卫星时间', '惯性系姿态角 X', '惯性系姿态角 Y', '惯性系姿态角 Z',
                                     '轨道系姿态角 X', '轨道系姿态角 Y', '轨道系姿态角 Z'],
                            '0111': ['rcv_time', 'eps_time']}

        # sep dictionary
        sep_dict = {'0811': ',', '0111': '\s+'}

        # skiprow
        skiprows_dict = {'0811': 2, '0111': 31}
        self.df = dd.read_csv(urlpath=self.path_list, sep=sep_dict[flag], header=None,
                              engine='c', skiprows=skiprows_dict[flag],
                              storage_options=dict(auto_mkdir=False), names=header_dict[flag],
                              dtype=dtype_dict[flag], encoding='gb2312')
        self.df = self.df[used_column_dict[flag]]
        self.df = self.df.rename(columns=dict(zip(self.df.columns, new_header_dict[flag])))
        self.df = self.df.drop_duplicates(subset=['rcv_time'])

        return self.df


def plot_freq_response(taps, fq):
    """
    Plot the magnitude response of the filter
    """
    plt.style.use(['science', 'no-latex', 'high-vis'])
    fig, ax = plt.subplots(figsize=(15, 8))
    plt.clf()
    w, h = freqz(taps, worN=8000)
    plt.loglog(w / np.pi * fq, abs(h), linewidth=2)
    plt.xlabel('Frequency [Hz]', fontsize=15)
    plt.ylabel('Gain', fontsize=15)
    plt.title('Frequency Response', fontsize=20)
    plt.tick_params(labelsize=25, width=2.9)
    plt.gca().spines['left'].set_linewidth(2)
    plt.gca().spines['top'].set_linewidth(2)
    plt.gca().spines['right'].set_linewidth(2)
    plt.gca().spines['bottom'].set_linewidth(2)
    plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    plt.show()


# check leap year
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


def rcvt2gpst(df_rcv, df_eps):
    """ This function is designed to transform satellite time (receiving time) to universal gps
        time

    Args:
        df_rcv (dask.dataframe): dataframe containing the rcv_time column
        df_eps (dask.dataframe): dataframe containing the time offset of rcv_time from gps_time
    """
    eps_time_interp = interpolate.interp1d(df_eps['rcv_time'].compute().to_numpy(),
                                           df_eps['eps_time'].compute().to_numpy(), kind='linear',
                                           bounds_error=False, fill_value='extrapolate')
    eps_time = eps_time_interp(df_rcv.self_rcv_time.compute().to_numpy())
    df_rcv = df_rcv.reset_index().set_index('index')
    df_rcv['gps_time'] = df_rcv.self_rcv_time.compute() + eps_time + 52588800.0

    # resort according to the gps_time
    df_rcv = df_rcv.nsmallest(df_rcv.__len__(), 'gps_time')

    return df_rcv


def get_run_time(func):
    def call_func(*args, **kwargs):
        begin_time = time.time()
        ret = func(*args, **kwargs)
        end_time = time.time()
        run_time = end_time - begin_time
        print(str(func.__name__) + "函数运行时间为" + str(run_time))
        return ret
    return call_func


# main function
@get_run_time
def main(year, month):
    """ This is a function to wrap the attitude data
    """

    null_df = dd.from_pandas(pd.DataFrame([]), npartitions=1)
    # read attitude files
    attd_file = InputFile(flag='0811', year=year, month=month, dir_path='..//input//',
                          df=null_df, path_list=[], output_filename='', days=0)
    attd_file.path_list, attd_file.output_filename, attd_file.days = attd_file.io_file_paths()
    attd_file.df = attd_file.read_data('0811')

    # read time offset files
    toff_file = InputFile(flag='0111', year=year, month=month, dir_path='..//input//',
                          df=null_df, path_list=[], output_filename='', days=0)
    toff_file.path_list, toff_file.output_filename, toff_file.days = toff_file.io_file_paths()
    toff_file.df = toff_file.read_data('0111')

    # creat an instance
    taiji_01_att = Att(attd_file.df, attd_file.days)
    taiji_01_att.date_reform('yyyy-mm-ddTHH:MM:SS')
    taiji_01_att.df = rcvt2gpst(taiji_01_att.df, toff_file.df)
    flags_list = taiji_01_att.equidistant_quantity(['igrf_eul_x', 'igrf_eul_y', 'igrf_eul_z',
                                                    'tf_eul_x', 'tf_eul_y', 'tf_eul_z'], 0.25)
    filtered_list = taiji_01_att.down_sample(flags_list, cutoff_hz=1., fq=0.5 / 0.25)
    # taiji_01_att.write_output_file_header(attd_file.output_filename)
    # taiji_01_att.write_file(filtered_list, attd_file.output_filename)

    # freq_igrf_x, psd_igrf_x = welch(
    #     taiji_01_att.df.igrf_eul_x.compute().to_numpy(), 4., 'hanning', taiji_01_att.df.__len__(), scaling='density')
    # freq_igrf_y, psd_igrf_y = welch(
    #     taiji_01_att.df.igrf_eul_y.compute().to_numpy(), 4., 'hanning', taiji_01_att.df.__len__(), scaling='density')
    # freq_igrf_z, psd_igrf_z = welch(
    #     taiji_01_att.df.igrf_eul_z.compute().to_numpy(), 4., 'hanning', taiji_01_att.df.__len__(), scaling='density')
# 
    # plt.style.use(['science', 'no-latex', 'high-vis'])
    # fig, ax = plt.subplots(figsize=(15, 8))
    # plt.plot(taiji_01_att.df.gps_time.compute().to_numpy(), taiji_01_att.df.igrf_eul_x.compute().to_numpy(),
    #          linewidth=2, label='gcrs2srf_x')
    # plt.plot(taiji_01_att.df.gps_time.compute().to_numpy(), taiji_01_att.df.igrf_eul_y.compute().to_numpy(),
    #          linewidth=2, label='gcrs2srf_y')
    # plt.plot(taiji_01_att.df.gps_time.compute().to_numpy(), taiji_01_att.df.igrf_eul_z.compute().to_numpy(),
    #          linewidth=2, label='gcrs2srf_z')
    # plt.tick_params(labelsize=25, width=2.9)
    # ax.yaxis.get_offset_text().set_fontsize(24)
    # ax.xaxis.get_offset_text().set_fontsize(24)
    # plt.xlabel('GPS time [s]', fontsize=20)
    # plt.ylabel('Attitude angle [degree]', fontsize=20)
    # plt.legend(fontsize=20, loc='best')
    # plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    # plt.gca().spines['left'].set_linewidth(2)
    # plt.gca().spines['top'].set_linewidth(2)
    # plt.gca().spines['right'].set_linewidth(2)
    # plt.gca().spines['bottom'].set_linewidth(2)
    # 
    # fig, ax = plt.subplots(figsize=(15, 8))
    # plt.loglog(freq_igrf_x, np.sqrt(psd_igrf_x),
    #          linewidth=2, label='gcrs2srf_x')
    # plt.loglog(freq_igrf_y, np.sqrt(psd_igrf_y),
    #          linewidth=2, label='gcrs2srf_y')
    # plt.loglog(freq_igrf_z, np.sqrt(psd_igrf_z),
    #          linewidth=2, label='gcrs2srf_z')
    # plt.tick_params(labelsize=25, width=2.9)
    # ax.yaxis.get_offset_text().set_fontsize(24)
    # ax.xaxis.get_offset_text().set_fontsize(24)
    # plt.xlabel('Frequency [Hz]', fontsize=20)
    # plt.ylabel(r'$PSD \quad [deg/\sqrt{hz}]$', fontsize=20)
    # plt.legend(fontsize=20, loc='best')
    # plt.grid(True, which='both', ls='dashed', color='0.5', linewidth=0.6)
    # plt.gca().spines['left'].set_linewidth(2)
    # plt.gca().spines['top'].set_linewidth(2)
    # plt.gca().spines['right'].set_linewidth(2)
    # plt.gca().spines['bottom'].set_linewidth(2)
    # 
    # plt.show()


if __name__ == '__main__':
    main(2019, 9)
