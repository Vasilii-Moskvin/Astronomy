import csv
import os
import pylab
import numpy as np
from collections import namedtuple, deque
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
#import seaborn as sns
import math
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 1000000000


LightCurve = namedtuple('LightCurve', ('time', 'ampl', 'err_ampl', 'smooth_ampl', 'smooth_err'))
LightCurve1 = namedtuple('LightCurve1', ('time', 'delta_diff_ampl'))

graph2D_titel = dict(DP='D’Agostino and Pearson test. p-value',
                     KS='Kolmogorov-Smirnov test. p-value',
                     SW='Shapiro-Wilk test. p-value',
                     NRW='Normalized residuals')


def open_csv(file_path):
    data = []
    with open(file_path, 'r') as f:
        header = f.readline().strip().split(',')

    with open(file_path, 'r') as csv_file:
        csv_file = csv.DictReader(csv_file, delimiter=',')
        raw_data = [row for row in csv_file]

    for src in raw_data:
        temp_dct = dict()
        for key, value in src.items():
            temp_dct[key] = float(value)
        data.append(temp_dct)

    data = list(filter(lambda x: x['delta'] != 0, data))

    return header, data


def column_from_csv(data, column_name):
    return list(map(lambda x: x[column_name], data))


def lc_data(data):
    time = column_from_csv(data, 'time')
    #ampl = list(sg.detrend(column_from_csv(data, 'ampl'), type='constant'))
    ampl = column_from_csv(data, 'delta')
    err_ampl = column_from_csv(data, 'err_delta')

    return time, ampl, err_ampl


def smooth_ampl_fun(ampl, err, n):
    ampl_1 = list(ampl)
    err_1 = list(err)
    smooth_ampl = ampl_1[:int(n)]
    smooth_err = err_1[:int(n)]
    temp_smooth_ampl = deque(ampl_1[:int(n)])
    temp_smooth_err = deque(err_1[:int(n)])
    ampl_local, sigma_local = mean_std(temp_smooth_ampl, temp_smooth_err)
    length = len(ampl)
    for i in range(int(n), length, 1):
        if i < length - 2:
            if ampl_1[i] + err_1[i] >= ampl_local + 2 * sigma_local:
                if ampl_1[i + 1] + err_1[i + 1] >= ampl_local + 2 * sigma_local and ampl_1[i + 2] + err_1[i + 2] >= ampl_local +  2 * sigma_local:
                    pass
                else:
                    temp_smooth_ampl.popleft()
                    temp_smooth_ampl.append(ampl_1[i])

                    temp_smooth_err.popleft()
                    temp_smooth_err.append(err_1[i])

                    ampl_local, sigma_local = mean_std(temp_smooth_ampl, temp_smooth_err)
            else:
                temp_smooth_ampl.popleft()
                temp_smooth_ampl.append(ampl_1[i])

                temp_smooth_err.popleft()
                temp_smooth_err.append(err_1[i])

                ampl_local, sigma_local = mean_std(temp_smooth_ampl, temp_smooth_err)
        else:
            ampl_local, sigma_local = mean_std(temp_smooth_ampl, temp_smooth_err)

        smooth_ampl.append(ampl_local)
        smooth_err.append(sigma_local)

    return np.array(smooth_ampl), np.array(smooth_err)


def split_smoth_lc(time, diff_smooth):
    time_flare_arr = []
    flare_arr = []
    time_check_arr = []
    check_arr = []
    sigma = np.std(diff_smooth)
    for src in zip(time, diff_smooth):
        if math.fabs(src[1]) <= 2 * sigma:
            time_check_arr.append(src[0])
            check_arr.append(src[1])
        else:
            time_flare_arr.append(src[0])
            flare_arr.append(src[1])

    return time_flare_arr, flare_arr, time_check_arr, check_arr


def mean_std0(x_i, err_i):
    x_i = np.array(x_i)
    err_i = np.array(err_i)

    n = len(x_i)
    # C = np.mean(err_i)

    # p_i = (C ** 2) / (err_i ** 2)
    p_i = 1 / err_i
    p = sum(p_i)

    mean_x0 = (1 / p) * sum(p_i * x_i)
    # sigma_eps = C / math.sqrt(p)
    # sigma_star = math.sqrt((1 / (n - 1)) * (sum(p_i * x_i ** 2)))# - p * mean_x0 ** 2))
    std_output = math.sqrt((1 / (n - 1)) * (sum(p_i * x_i ** 2)))
    # sigma_star_eps = sigma_star / math.sqrt(p)

    # if sigma_eps > sigma_star_eps:
    #     std_output = 0.5 * (sigma_eps + sigma_star_eps)
    # else:
    #     std_output = sigma_star_eps

    std_x0 = std_output / math.sqrt(n)

    return mean_x0, std_x0, std_output


def mean_std(x_i, err_i):
    mean_x0, std_x0, std_output = mean_std0(x_i, err_i)

    return mean_x0, std_x0 + std_output


def draw_LSNR(time, ampl, err_ampl, smooth_ampl, smooth_err, delta_diff_ampl, flare_list):
    number = 2
    f, axarr = plt.subplots(number, sharex=True)
    axarr[0].plot(time, ampl, linestyle=' ', marker='.', color='red', markersize=1)
    axarr[0].fill_between(time, ampl + err_ampl, ampl - err_ampl, color='grey', alpha=0.3)
    axarr[0].plot(time, smooth_ampl, linestyle=' ', marker='.', color='blue', markersize=1)
    axarr[0].fill_between(time,  smooth_ampl + 3 * smooth_err, smooth_ampl - 3 * smooth_err, color='pink', alpha=0.3)
    axarr[0].fill_between(time,  smooth_ampl + 2 * smooth_err, smooth_ampl - 2 * smooth_err, color='green', alpha=0.3)
    axarr[0].fill_between(time, smooth_ampl + smooth_err, smooth_ampl - smooth_err, color='blue', alpha=0.3)
    for src in flare_list:
        axarr[0].axvspan(src[0].time, src[1].time, alpha=0.3, color='blue')
    axarr[0].set_title('Light and smooth curve')

    axarr[1].plot(time, delta_diff_ampl, linestyle=' ', marker='.', color='red', markersize=1)
    axarr[1].set_title('Normalized residuals', y=0.99)
    axarr[1].set_xlabel('Time, sec')

    f.text(0.03, 0.5, 'Intensity, count/sec', va='center', rotation='vertical')

    for index in range(number):
        if index != 2 and index != 3:
            axarr[index].grid(color='grey', linestyle='--', linewidth=0.5)

    savename = 'C:\\Users\\vasil\\Desktop\\{}'.format('LSNR')
    plt.savefig(savename,
                  dpi=500,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None,
                  fmt='png')

    plt.close()


def draw_four_2D_graph(df, df_norm, df_kstest, df_shapiro):
    #number = 4
    f, axarr = plt.subplots(2, 2, sharex='col', sharey='row')

    xmin = df.columns[0]
    xmax = df.columns[-1]
    ymin = df.index[0]
    ymax = df.index[-1]
    z00_plot = axarr[0][0].matshow(df, aspect='auto', extent=(xmin, xmax, ymax, ymin))
    axarr[0][0].xaxis.set_label_position('bottom')
    axarr[0][0].xaxis.set_ticks_position('bottom')
    axarr[0][0].invert_yaxis()
    plt.colorbar(z00_plot, ax=axarr[0][0])
    axarr[0][0].set_title('Normalized residuals')

    z01_plot = axarr[0][1].matshow(df_norm, aspect='auto', extent=(xmin, xmax, ymax, ymin))
    axarr[0][1].xaxis.set_label_position('bottom')
    axarr[0][1].xaxis.set_ticks_position('bottom')
    axarr[0][1].invert_yaxis()
    plt.colorbar(z01_plot, ax=axarr[0][1])
    axarr[0][1].set_title('D’Agostino and Pearson. p-value')

    z10_plot = axarr[1][0].matshow(df_kstest, aspect='auto', extent=(xmin, xmax, ymax, ymin))
    axarr[1][0].xaxis.set_label_position('bottom')
    axarr[1][0].xaxis.set_ticks_position('bottom')
    axarr[1][0].invert_yaxis()
    axarr[1][0].set_xlabel('Time, sec')
    plt.colorbar(z10_plot, ax=axarr[1][0])
    axarr[1][0].set_title('Kolmogorov-Smirnov. p-value', y=0.99)

    z11_plot = axarr[1][1].matshow(df_shapiro, aspect='auto', extent=(xmin, xmax, ymax, ymin))
    axarr[1][1].xaxis.set_label_position('bottom')
    axarr[1][1].xaxis.set_ticks_position('bottom')
    axarr[1][1].invert_yaxis()
    axarr[1][1].set_xlabel('Time, sec')
    plt.colorbar(z11_plot, ax=axarr[1][1])
    axarr[1][1].set_title('Shapiro-Wilk. p-value', y=0.99)

    # f.text(0.5, 0.48, 'p-value within widows in normalized residuals', ha='center')
    # f.text(0.5, 0.02, 'Time, sec', ha='center')
    f.text(0.03, 0.5, 'Windows size, sec', va='center', rotation='vertical')

    savename = 'C:\\Users\\vasil\\Desktop\\{}'.format('NDKS')
    plt.savefig(savename,
                dpi=500,
                facecolor='w',
                edgecolor='w',
                orientation='portrait',
                parpertype=None,
                format=None,
                transparent=False,
                bbox_inches=None,
                pad_inches=0.1,
                frameon=None,
                fmt='png')

    #plt.show()
    plt.close()


def draw_two_2D_graph(df0, df1, graph_type0, graph_type1):
    number = 2
    f, axarr = plt.subplots(number, sharex=True)

    xmin = df0.columns[0]
    xmax = df0.columns[-1]
    ymin = df0.index[0]
    ymax = df0.index[-1]
    z0_plot = axarr[0].matshow(df0, aspect='auto', extent=(xmin, xmax, ymax, ymin))
    axarr[0].xaxis.set_label_position('bottom')
    axarr[0].xaxis.set_ticks_position('bottom')
    axarr[0].invert_yaxis()
    cbar0_ax = f.add_axes([0.92, 0.53, 0.01, 0.35])
    f.colorbar(z0_plot, cax=cbar0_ax)
    # plt.colorbar(z0_plot, ax=axarr[0])
    axarr[0].set_title('{}'.format(graph2D_titel[graph_type0]))

    xmin = df1.columns[0]
    xmax = df1.columns[-1]
    ymin = df1.index[0]
    ymax = df1.index[-1]
    z1_plot = axarr[1].matshow(df1, aspect='auto', extent=(xmin, xmax, ymax, ymin))
    axarr[1].xaxis.set_label_position('bottom')
    axarr[1].xaxis.set_ticks_position('bottom')
    axarr[1].invert_yaxis()
    axarr[1].set_xlabel('Time, sec')
    cbar1_ax = f.add_axes([0.92, 0.11, 0.01, 0.35])
    f.colorbar(z1_plot, cax=cbar1_ax)
    axarr[1].set_title('{}'.format(graph2D_titel[graph_type1]), y=0.99)

    # f.text(0.5, 0.48, 'p-value within widows in normalized residuals', ha='center')
    # f.text(0.5, 0.02, 'Time, sec', ha='center')
    f.text(0.03, 0.5, 'Windows size, sec', va='center', rotation='vertical')

    savename = 'C:\\Users\\vasil\\Desktop\\{}_{}'.format(graph_type0, graph_type1)
    plt.savefig(savename,
                dpi=500,
                facecolor='w',
                edgecolor='w',
                orientation='portrait',
                parpertype=None,
                format=None,
                transparent=False,
                bbox_inches=None,
                pad_inches=0.1,
                frameon=None,
                fmt='png')

    # plt.show()
    plt.close()


def draw_one_2D_graph(df, graph_type):

    xmin = df.columns[0]
    xmax = df.columns[-1]
    ymin = df.index[0]
    ymax = df.index[-1]

    plt.matshow(df, aspect='auto', extent=(xmin, xmax, ymax, ymin))
    plt.colorbar()

    plt.gca().xaxis.set_label_position('bottom')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().invert_yaxis()
    plt.gca().set_title('{}'.format(graph2D_titel[graph_type]))
    plt.gca().set_xlabel('Time, sec')
    plt.gca().set_ylabel('Windows size, sec')


    savename = 'C:\\Users\\vasil\\Desktop\\{}'.format(graph_type)
    plt.savefig(savename,
                dpi=500,
                facecolor='w',
                edgecolor='w',
                orientation='portrait',
                parpertype=None,
                format=None,
                transparent=False,
                bbox_inches=None,
                pad_inches=0.1,
                frameon=None,
                fmt='png')

    #plt.show()
    plt.close()


def draw_one_1Dgraph(x, y, err_y, graph_type):
    plt.plot(x, y, linestyle=' ', marker='.', color='red', markersize=1)
    if graph_type != 'NR':
        plt.fill_between(x, y + err_y, y - err_y, color='grey', alpha=0.3)
        if graph_type == 'LC':
            plt.gca().set_title('Light curve')
        else:
            plt.gca().set_title('Smooth light curve')
    else:
        plt.gca().set_title('Normalized residuals')
    plt.gca().set_xlabel('Time, sec')
    plt.gca().set_ylabel('Intensity, count/sec')
    plt.gca().grid(color='grey', linestyle='--', linewidth=0.5)

    savename = 'C:\\Users\\vasil\\Desktop\\{}'.format(graph_type)
    plt.savefig(savename,
                dpi=500,
                facecolor='w',
                edgecolor='w',
                orientation='portrait',
                parpertype=None,
                format=None,
                transparent=False,
                bbox_inches=None,
                pad_inches=0.1,
                frameon=None,
                fmt='png')

    #plt.show()
    plt.close()


def draw_LSC(x, y, err_y, ys, err_ys):
    plt.plot(x, y, linestyle=' ', marker='.', color='red', markersize=1)
    plt.fill_between(x, y + err_y, y - err_y, color='grey', alpha=0.3)
    plt.plot(x, ys, linestyle=' ', marker='.', color='blue', markersize=1)
    plt.fill_between(x, ys + err_ys, ys - err_ys, color='blue', alpha=0.3)

    plt.gca().set_title('Light and smooth curves')
    plt.gca().set_xlabel('Time, sec')
    plt.gca().set_ylabel('Intensity, count/sec')
    plt.gca().grid(color='grey', linestyle='--', linewidth=0.5)

    savename = 'C:\\Users\\vasil\\Desktop\\{}'.format('LSC')
    plt.savefig(savename,
                dpi=500,
                facecolor='w',
                edgecolor='w',
                orientation='portrait',
                parpertype=None,
                format=None,
                transparent=False,
                bbox_inches=None,
                pad_inches=0.1,
                frameon=None,
                fmt='png')

    #plt.show()
    plt.close()
    

def fun(time, tau_1, tau_2, ts, height):
    time_new = list(filter(lambda x: x > ts, time))
    lam = math.sqrt(math.exp(2 * (tau_1 / tau_2)))
    i_output = list(map(lambda t: height * lam * math.exp(-((tau_1 / (t - ts)) + ((t - ts) / tau_2))), time_new))

    return time_new, i_output


def diff_ampl_fun(ampl, smooth_ampl, err_ampl, smooth_err):
    diff_ampl = ampl - smooth_ampl
    diff_err = err_ampl + smooth_err
    delta_diff_ampl = diff_ampl / diff_err
    
    return diff_ampl, diff_err, delta_diff_ampl


def normaltest_fun(delta_diff_ampl, n):
    nt = np.zeros(len(delta_diff_ampl[:int(n / 2)]))
    for index, value in enumerate(delta_diff_ampl):
        if index >= int(n / 2) and index < len(delta_diff_ampl) - int(n / 2):
            temp_nt = delta_diff_ampl[index - int(n / 2):index + int(n / 2): 1]
            new_nt = stats.normaltest(temp_nt).pvalue
            nt = np.append(nt, new_nt)

    nt = np.append(nt, np.zeros(len(delta_diff_ampl[int(-n / 2):])))
    
    return nt


def kstest_fun(delta_diff_ampl, n):
    nt = np.zeros(len(delta_diff_ampl[:int(n / 2)]))
    for index, value in enumerate(delta_diff_ampl):
        if index >= int(n / 2) and index < len(delta_diff_ampl) - int(n / 2):
            temp_nt = delta_diff_ampl[index - int(n / 2):index + int(n / 2): 1]
            new_nt = stats.kstest(temp_nt, 'norm').pvalue
            nt = np.append(nt, new_nt)

    nt = np.append(nt, np.zeros(len(delta_diff_ampl[int(-n / 2):])))

    return nt


def shapiro_fun(delta_diff_ampl, n):
    nt = np.zeros(len(delta_diff_ampl[:int(n / 2)]))
    for index, value in enumerate(delta_diff_ampl):
        if index >= int(n / 2) and index < len(delta_diff_ampl) - int(n / 2):
            temp_nt = delta_diff_ampl[index - int(n / 2):index + int(n / 2): 1]
            new_nt = stats.shapiro(temp_nt)[1]
            nt = np.append(nt, new_nt)

    nt = np.append(nt, np.zeros(len(delta_diff_ampl[int(-n / 2):])))

    return nt


def check_is_flare(tmp):
    sigma2 = deque()
    for src in tmp:
        # print('{}: {}'.format(src.time, (src.ampl + src.err_ampl) / (src.smooth_ampl + src.smooth_err) ))
        if src.ampl + src.err_ampl >=  src.smooth_ampl + 1.5 * src.smooth_err:
            sigma2.append(src)
            if len(sigma2) > 2:
                if filter(lambda x: x.ampl + x.err_ampl >= x.smooth_ampl + 2 * x.smooth_err, sigma2):
                    return True
        else:
            if sigma2:
                sigma2 = deque()
    return False


def check_is_flare1(tmp):
    sigma2 = deque()
    for src in tmp:
        # print('{}: {}'.format(src.time, (src.ampl + src.err_ampl) / (src.smooth_ampl + src.smooth_err) ))
        if src.delta_diff_ampl >=  1.5:
            sigma2.append(src)
            if len(sigma2) > 2:
                if filter(lambda x: x.delta_diff_ampl >= 2, sigma2):
                    return True
        # else:
        #     if sigma2:
        #         sigma2 = deque()
    return False


def get_flares(time, ampl, err_ampl, smooth_ampl, smooth_err):
    flare_list = []
    xren = zip(time, ampl, err_ampl, smooth_ampl, smooth_err)
    xren_lst = [LightCurve(*src) for src in xren]
    tmp = []
    rec_flare = False
    deq = deque()
    flag = False
    for src in xren_lst:
        if src.ampl + src.err_ampl > src.smooth_ampl: #+ 0.5 * src.smooth_err:
            if not rec_flare:
                rec_flare = True
                deq.append(src)
        else:
            if rec_flare:
                deq.append(src)
                rec_flare = False
                tmp.append(src)
                is_flare = check_is_flare(tmp)
                tmp = []
                if is_flare:
                    flare_list.append(tuple(deq))
                deq = deque()
        if rec_flare:
            tmp.append(src)

    for src in flare_list:
        print('{}\n{}'.format(*src))


def get_flares1(time, delta_diff_ampl):
    flare_list = []
    xren = zip(time, delta_diff_ampl)
    xren_lst = [LightCurve1(*src) for src in xren]
    tmp = []
    rec_flare = False
    deq = deque()
    flag = False
    for src in xren_lst:
        if src.delta_diff_ampl > 0:
            if not rec_flare:
                rec_flare = True
                deq.append(src)
        else:
            if rec_flare:
                deq.append(src)
                rec_flare = False
                tmp.append(src)
                is_flare = check_is_flare1(tmp)
                tmp = []
                if is_flare:
                    flare_list.append(tuple(deq))
                deq = deque()
        if rec_flare:
            tmp.append(src)

    for src in flare_list:
        print('{}\n{}'.format(*src))

    return flare_list


def calc_sigma_noise(time, ampl, time_start, time_end):
	'''
		Функция расчитывает среднее значение и стандартное отклонение шума 
		на выбранном интервале [time_start, time_end]
	'''
	time_ampl_all = zip(time, ampl)
    time_ampl_in_interval = [(t, a) for t, a in time_ampl_all if time_start <= t <= time_end]
    ampl_in_interval = list(map(lambda x: x[1], time_ampl_in_interval))
    average_noise_intensity, std_noise_intensity = calc_average_std_intensity(ampl_in_interval)


	return average_noise_intensity, std_noise_intensity

def calc_average_std_intensity(ampl):
	'''
		Функция возвращает среднее значение и стандартное отклонение величины ampl
	'''
	n = len(ampl)
    average_ampl = sum(ampl) / n

    element_of_sum_for_std = [(a - average_ampl) ** 2 for a in ampl]
    std_ampl = math.sqrt(sum(element_of_sum_for_std) / (n - 1))

    return average_ampl, std_ampl 



def main():
    #lc_file_path = os.path.abspath(input('Enter path to file with light curve data:\n'))
    lc_file_path = os.path.abspath(r'C:\Users\vasil\Desktop\results_2\CN_Leo\090128_04\090128_04_1sec_lc.csv')
    header, data = open_csv(lc_file_path)
    time, ampl, err_ampl = map(np.array, lc_data(data))
    # df = pd.DataFrame(columns=time)
    # df_DP = pd.DataFrame(columns=time)
    # df_KS = pd.DataFrame(columns=time)
    # df_SW = pd.DataFrame(columns=time)
    
    n = 600
    smooth_ampl, smooth_err = smooth_ampl_fun(ampl, err_ampl, n)

    diff_ampl, diff_err, delta_diff_ampl = diff_ampl_fun(ampl, smooth_ampl, smooth_err, err_ampl)
    flare_list = get_flares1(time, delta_diff_ampl)#, err_ampl, smooth_ampl, smooth_err)
    # df.loc[n] = delta_diff_ampl
        # df_DP.loc[n] = normaltest_fun(delta_diff_ampl, n)
        # df_KS.loc[n] = kstest_fun(delta_diff_ampl, n)
        # df_SW.loc[n] = shapiro_fun(delta_diff_ampl, n)
    print(n)
        #print('normal skewtest teststat = {} pvalue = {}'.format(*stats.normaltest(delta_diff_ampl)))
        # if n == 40:
    smooth_ampl0 = smooth_ampl[:]
    smooth_err0 = smooth_err[:]
    diff_ampl0 = diff_ampl[:]
    diff_err0 = diff_err[:]
    delta_diff_ampl0 = delta_diff_ampl[:]
    # for n in range(39, 41):
    #     #n = 40
    #     smooth_ampl, smooth_err = smooth_ampl_fun(ampl, err_ampl, n)
    #     diff_ampl, diff_err, delta_diff_ampl = diff_ampl_fun(ampl, smooth_ampl, smooth_err, err_ampl)
    #     df.loc[n] = delta_diff_ampl
    #     # df_DP.loc[n] = normaltest_fun(delta_diff_ampl, n)
    #     # df_KS.loc[n] = kstest_fun(delta_diff_ampl, n)
    #     # df_SW.loc[n] = shapiro_fun(delta_diff_ampl, n)
    #     print(n)
    #     #print('normal skewtest teststat = {} pvalue = {}'.format(*stats.normaltest(delta_diff_ampl)))
    #     if n == 40:
    #         smooth_ampl0 = smooth_ampl[:]
    #         smooth_err0 = smooth_err[:]
    #         diff_ampl0 = diff_ampl[:]
    #         diff_err0 = diff_err[:]
    #         delta_diff_ampl0 = delta_diff_ampl[:]


    draw_LSNR(time, ampl, err_ampl, smooth_ampl0, smooth_err0, delta_diff_ampl0, flare_list)
    draw_one_1Dgraph(time, ampl, err_ampl, 'LC')
    draw_one_1Dgraph(time, smooth_ampl0, smooth_err0, 'SC')
    draw_one_1Dgraph(time, delta_diff_ampl0, diff_err0, 'NR')
    draw_LSC(time, ampl, err_ampl, smooth_ampl0, smooth_err0)
    # draw_one_2D_graph(df, 'NRW')
    # draw_one_2D_graph(df_DP, 'DP')
    # draw_one_2D_graph(df_KS, 'KS')
    # draw_one_2D_graph(df_SW, 'SW')
    # draw_two_2D_graph(df, df_DP, 'NRW', 'DP')
    # draw_two_2D_graph(df, df_KS, 'NRW', 'KS')
    # draw_two_2D_graph(df, df_SW, 'NRW', 'SW')
    # draw_four_2D_graph(df, df_DP, df_KS, df_SW)

    #print(df.index)
    #sns.heatmap(df, square=True, linewidths=.5, cbar_kws={"shrink": .5})


if __name__ == '__main__':
    main()