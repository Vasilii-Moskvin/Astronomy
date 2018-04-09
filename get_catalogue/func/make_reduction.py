import copy
import os
import numpy as np
from func.StarFromCSV import StarFromCSV
from func.work_with_csv import get_column
from collections import namedtuple
from func.graphics import save_regress


ABC = namedtuple('ABC', ('a', 'b', 'c'))


def get_ABC(to_sys_data_corr):
    '''
    Gets the coefficients obtained by the method of least squares.
    :param to_sys_data_corr: processed data from the compilation catalogue.
    :return: the coefficients obtained by the method of least squares.
    '''
    ABC_by_filt = dict(map(lambda key: (key, []), StarFromCSV.filts))

    if len(StarFromCSV.filts) == 1:
        filt = StarFromCSV.filts[0]

        best_data = list(filter(lambda i: i['{}_std_mag'.format(filt)] < 0.03, to_sys_data_corr))
        x = get_column(best_data, '{}_avr_mag'.format(filt))
        y = get_column(best_data, '{}_mag_cat'.format(filt))

        A = np.vstack([x, np.ones(len(x))]).T
        k, b = np.linalg.lstsq(A, y)[0]

        ABC_by_filt[filt] = ABC(k, 0, b)
    else:
        for filt in StarFromCSV.filts:
            best_data = list(filter(lambda i: i['{}_std_mag'.format(filt)] < 0.03, to_sys_data_corr))

            filt2 = StarFromCSV.choose_filter(filt)
            f1_inst, f2_f1_inst, f1_cat = div_by_filter(best_data, filt, filt2)

            A = np.vstack([f1_inst, f2_f1_inst, np.ones(len(f1_inst))]).T
            a, b, c = np.linalg.lstsq(A, f1_cat)[0]
            ABC_by_filt[filt] = ABC(a, b, c)

    return ABC_by_filt


def div_by_filter(data, filt1, filt2):
    '''
    Collects information about the value of stellar magnitudes from the processed data from the compilation catalogue.
    The values are returned: the instrumental stellar magnitudes in a certain filter,
    the difference between these magnitudes and the corresponding stellar magnitudes in another filter,
    the stellar magnitudes in a certain filter in the reference catalogue.
    :param data: the processed data from the compilation catalogue.
    :param filt1: the selected filter.
    :param filt2: filter, the difference in stellar magnitudes of which and star values in the selected filter 
    will be determined by the color index.
    :return: the instrumental stellar magnitudes in a certain filter,
    the difference between these magnitudes and the corresponding stellar magnitudes in another filter,
    the stellar magnitudes in a certain filter in the reference catalogue.
    '''
    filt1_inst = []
    f2_f1_inst = []
    filt1_cat = []

    for src in data:
        filt1_inst.append(src['{}_avr_mag'.format(filt1)])
        f2_f1_inst.append(src['{}_avr_mag'.format(filt2)] - src['{}_avr_mag'.format(filt1)])
        filt1_cat.append(src['{}_mag_cat'.format(filt1)])

    return filt1_inst, f2_f1_inst, filt1_cat


def make_reducted_catalogue(to_sys_data_corr):
    '''
    Makes reducted catalogue.
    :param ABC_by_filt: the coefficients obtained by the method of least squares, depending on the filters.
    :param to_sys_data_corr: catalogue with information about satrs.
    :return: the reducted catalogue.
    '''
    ABC_by_filt = get_ABC(to_sys_data_corr)
    make_data = copy.deepcopy(to_sys_data_corr)

    if len(StarFromCSV.filts) == 1:
        for star in make_data:
            star.reduction_AC(ABC_by_filt)
    else:
        for star in make_data:
            star.reduction_ABC(ABC_by_filt)

    return make_data, ABC_by_filt


def make_reducted_catalogue_1(data, ABC_by_filt, two=False):
    if two:
        make_data = list(filter(lambda star: all(list(map(lambda x: int(star[x]) > 2, StarFromCSV.n_fields))), copy.deepcopy(data)))
    else:
        make_data = copy.deepcopy(data)

    if len(StarFromCSV.filts) == 1:
        for star in make_data:
            star.reduction_AC(ABC_by_filt)
    else:
        # for star in make_data:
        #     star.reduction_ABC(ABC_by_filt)
        return None

    return make_data


def save_regress_line(dir_path, reducted_catalogue, ABC_by_filt):
    '''
    Prepareses data and settings for saving regression lines. Saves regression lines.
    :param dir_path: path to the work folder.
    :param reducted_catalogue: the reducted catalogue.
    :param ABC_by_filt: the coefficients obtained by the method of least squares, depending on the filters.
    :return: saved images of regression lines.
    '''
    for filt in StarFromCSV.filts:
        x = get_column(reducted_catalogue, '{}_avr_mag'.format(filt))
        y = get_column(reducted_catalogue, '{}_mag_cat'.format(filt))

        A = np.vstack([x, np.ones(len(x))]).T
        k, b = np.linalg.lstsq(A, y)[0]

        y_mnk = [k * i + b for i in x]

        save_path = os.path.join(dir_path, 'regress_{}'.format(filt))
        x_lbl = 'Instrumental catalogue, mag'
        y_lbl = 'Catalogue, mag'
        tlt = 'Catalohue vs instrumental magnitudes ({} filter)'.format(filt)
        if ABC_by_filt[filt].b != 0:
            filt2 = StarFromCSV.choose_filter(filt)
            legend_lbl = '${}_c={:.3}*{}_i{:+.3}*({}_i-{}_i){:+.3}$\n$y={:.3}*x{:+.3}$'.format(filt, ABC_by_filt[filt].a, filt, ABC_by_filt[filt].b, filt2, filt, ABC_by_filt[filt].c, k, b)
        else:
            legend_lbl = '${}_c={:.3}*{}_i{:+.3}$\n$y={:.3}*x{:+.3}$'.format(filt, ABC_by_filt[filt].a, filt, ABC_by_filt[filt].c, k, b)
        save_regress(save_path, x, y, y_mnk, tlt, x_lbl, y_lbl, legend_lbl)
