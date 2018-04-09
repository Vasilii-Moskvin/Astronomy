import copy
import os
import numpy as np
from func.StarFromCSV import StarFromCSV
from pandas import DataFrame
from func.graphics import save_correlogram


def produce_corr_analysis(dir_path, data, to_sys_data, obj):
    '''
    Produces a correlation analysis of series of stellar magnitudes. Removes trends in the series of stellar magnitudes.
    :param dir_path: path to the work folder.
    :param data: data from the compilation catalogue.
    :param to_sys_data: processed data from the compilation catalogue.
    :param obj: coordinates of the observed object. This object will not be taken as a comparison object.
    :return: data with corrected trend. The first one for data from the compilation catalogue. 
    The second one for processed data from the compilation catalogue.
    '''
    #to_sys_data = smoothing(to_sys_data, 3)

    for filt in StarFromCSV.filts:
        ten_good = ten_best_star(dir_path, to_sys_data, obj, filt)
        stars_corr = stars_with_corr_coeff(ten_good)

        err_arr = get_err_arr(dir_path, to_sys_data, filt, stars_corr, ten_good)

        changedata(to_sys_data, err_arr, filt)
        changedata(data, err_arr, filt)

    return data, to_sys_data


def stars_with_corr_coeff(ten_good):
    '''
    Returns the list of stars with the correlation coefficient. 
    :param ten_good: the series of stellar magnitudes for 10 stars with a minimum error. 
    :return: the list of tuples in each of which located coordinates two stars and their correlation coefficient.
    '''
    rho = DataFrame(ten_good).corr()
    stars_corr = []
    for coord_1 in rho:
        for coord_2 in rho:
            if coord_1 != coord_2:
                if stars_corr:
                    for src in stars_corr:
                        if src[0] == coord_2 and src[1] == coord_1:
                            break
                        elif src[0] == coord_1 and src[1] == coord_2:
                            break
                    else:
                        stars_corr.append((coord_1, coord_2, rho[coord_1][coord_2]))
                else:
                    stars_corr.append((coord_1, coord_2, rho[coord_1][coord_2]))
    stars_corr.sort(key=lambda x: x[2], reverse=True)

    return stars_corr



def ten_best_star(dir_path, to_sys_data, obj, filt):
    '''
    Returns the series of stellar magnitudes (depending on the filter) for 10 stars with a minimum error. 
    Saves the resulting correlogram.
    :param dir_path: path to the folder where to save the correlogram.
    :param to_sys_data: processed data from the compilation catalogue.
    :param obj: coordinates of the observed object. This object will not be taken as a comparison object.
    :param filt: a filter in which 10 stars are looked for with minimal errors.
    :return: the series of stellar magnitudes (depending on the filter) for 10 stars with a minimum error. 
    The index is the coordinates of the objects.
    '''
    ten_stars = sorted(filter(lambda x: '{}_{}'.format(x['RAJ2000'], x['DEJ2000']) != obj, to_sys_data),
                          key=lambda x: x['{}_std_mag'.format(filt)])[:10]

    ten_good = dict()
    tensns = dict()
    for index, star in enumerate(ten_stars):
        temp = list(star.mag_digit_values(filt))
        ten_good['{}_{}'.format(star['RAJ2000'], star['DEJ2000'])] = temp
        tensns[index] = temp

    rhosns = DataFrame(tensns).corr()
    save_correlogram(dir_path, rhosns, filt)

    return ten_good


def smoothing(to_sys_data, n=3):
    '''
    Performs smoothing of data.
    :param to_sys_data: processed data from the compilation catalogue.
    :param n: Number of points using to smoothing (by default is 3).
    :return: smoothed of data.
    '''
    smth_to_sys_data = copy.deepcopy(to_sys_data)
    for star in smth_to_sys_data:
        for filt in star.mag_values.keys():
            new_values, err_new_vales = smoothing_run(list(star.mag_values[filt].values()), n)
            new_values_dct = dict(zip(star.mag_values[filt].keys(), new_values))
            for field in star.mag_values[filt].keys():
                star[field] = new_values_dct[field]

        star.update_mag()

    return smth_to_sys_data


def smoothing_run(lst, n=3):
    '''
    Starts smoothing of data.
    :param lst: list of values that need to be made smoothed.
    :param n: Number of points using to smoothing (by default is 3).
    :return: smoothed of data.
    '''
    new_lst = []
    error_lst = []
    for index, src in enumerate(lst):
        if index < len(lst) - n:
            temp = []
            for i in lst[index:index + n]:
                temp.append(i)
            new_lst.append(np.mean(temp))
            error_lst.append(np.std(temp))
        else:
            new_lst.append(src)
            error_lst.append(0)

    temp_lst = copy.deepcopy(new_lst[-1::-1])
    new_lst = []
    error_lst = []

    for index, src in enumerate(temp_lst):
        if index < len(temp_lst) - n:
            temp = []
            for i in temp_lst[index:index + n]:
                temp.append(i)
            new_lst.append(np.mean(temp))
            error_lst.append(np.std(temp))
        else:
            new_lst.append(src)
            error_lst.append(0)

    return new_lst[-1::-1], error_lst[-1::-1]


def get_err_arr(dir_path, to_sys_data, filt, stars_corr, ten_good):
    '''
    Returns the data for changes of the serieses of stellar magnitudes.
    :param dir_path: path to the work folder.
    :param to_sys_data: processed data from the compilation catalogue.
    :param filt: filter for which series of magnitude are made.
    :param stars_corr: the list of tuples in each of which located coordinates two stars and their correlation coefficient.
    :param ten_good: the series of stellar magnitudes for 10 stars with a minimum error.
    :return: the data for changes of the serieses of stellar magnitudes.
    '''
    err_arr = [0 for _ in StarFromCSV.mag_fields_name[filt]]
    savedata = []
    ref1_RADE = ''
    ref2_RADE = ''

    mean_old, std_old = get_mean_std(to_sys_data, filt)

    print('\n{}\nmean_old| std_old\n {} |{}'.format(filt, round(mean_old, 6), round(std_old, 6)))
    print('RA1        DE1       | RA2        DE2       | mean_std| std_std | corr')
    savedata.append('\n{}\nmean_old| std_old\n {} {}'.format(filt, round(mean_old, 6), round(std_old, 6)))
    #savedata.append('{} {} {}'.format(filt, round(mean_old, 6), round(std_old, 6)))
    savedata.append('RA1        DE1       | RA2        DE2       | mean_std| std_std | corr')

    for mkey1, mkey2, maxr in sorted(stars_corr, key=lambda x: x[2], reverse=True)[:5]:
        avr = np.mean(ten_good[mkey1])
        err_arr_0 = list(map(lambda x: x - avr, ten_good[mkey1]))

        temp_data = copy.deepcopy(to_sys_data)

        changedata(temp_data, err_arr_0, filt)
        mean_temp, std_temp = get_mean_std(temp_data, filt)

        if mean_temp + std_temp < mean_old + std_old:
            mean_old = mean_temp
            std_old = std_temp
            err_arr = err_arr_0
            savedata.append('{} {}| {} {}| {:.6f}| {:.6f}| {:.6f}| better'.format(*mkey1.split('_'),
                                                                                  *mkey2.split('_'),
                                                                                  mean_temp,
                                                                                  std_temp,
                                                                                  maxr))
            print('{} {}| {} {}| {:.6f}| {:.6f}| {:.6f}| better'.format(*mkey1.split('_'),
                                                                        *mkey2.split('_'),
                                                                        mean_temp,
                                                                        std_temp,
                                                                        maxr))
            ref1_RADE = mkey1
            ref2_RADE = mkey2
        else:
            savedata.append('{} {}| {} {}| {:.6f}| {:.6f}| {:.6f}'.format(*mkey1.split('_'),
                                                                          *mkey2.split('_'),
                                                                          mean_temp,
                                                                          std_temp,
                                                                          maxr))
            print('{} {}| {} {}| {:.6f}| {:.6f}| {:.6f}'.format(*mkey1.split('_'),
                                                                *mkey2.split('_'),
                                                                mean_temp,
                                                                std_temp,
                                                                maxr))
    with open(os.path.join(dir_path, 'log_{}.txt'.format(filt)), 'w', newline='') as f:
        for src in savedata:
            f.write('{}\n'.format(src))

    return err_arr


def changedata(data, err_arr, filt):
    '''
    Produces change_mag(err_arr, filt) over each star in the list data.
    :param data: the list of stars.
    :param err_arr: the data for changes.
    :param filt: filter for which series of magnitude are made.
    :return: change_mag(err_arr, filt) over each star in the list data.
    '''
    for star in data:
        star.change_mag(err_arr, filt)


def get_mean_std(data, filt):
    '''
    Calculates mean and std values for errors of all stars by the filter filt.
    :param data: the list of stars.
    :param filt: the error value for the magnitudes is selected for this filter.
    :return: mean and std vlaues for errors of all stars in the filter filt.
    '''
    stat = [src['{}_std_mag'.format(filt)] for src in data]
    mean = np.mean(stat)
    std = np.std(stat)
    return mean, std
