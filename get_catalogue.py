import os.path
import re
#import math
import numpy as np
import pyfits
from pandas import DataFrame
import csv
from collections import OrderedDict, namedtuple
from functools import reduce
import operator
import copy
import matplotlib
import matplotlib.pyplot as plt
from seaborn import heatmap
import json

ABC = namedtuple('ABC', ('a', 'b', 'c'))

class EquilPlaneError(BaseException):
    pass


class StarFromCSV(OrderedDict):
    header = []
    filts = []
    n_fields = OrderedDict()
    mag_fields_name = OrderedDict()
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.get_static_info(*args)
    
    @classmethod
    def reset_class_data(cls):
        cls.header = []
        cls.filts = []
        cls.n_fields = OrderedDict()
        cls.mag_fields_name = OrderedDict()

    @classmethod
    def get_static_info(cls, *args):
        if not cls.header:
            cls.header = [item[0] for item in args[0]]
        if not cls.filts:
            cls.filts = list(map(lambda x: re.findall(r'(^\w+)_n$', x)[0],
                              filter(lambda x: re.findall(r'^\w+_n$', x), cls.header)))
        if not cls.n_fields:
            cls.n_fields = OrderedDict(map(lambda field: (field, ''),
                                           filter(lambda x: re.findall(r'^\w+_n$', x), cls.header)))
        if not cls.mag_fields_name:
            for filt in cls.filts:
                cls.mag_fields_name[filt] = list(filter(lambda x: re.findall(r'{}_mag_\d+$'.format(filt), x),
                                                        cls.header))

    @classmethod
    def fill_n_fields(cls, dct_max_n_value):
        for field in cls.n_fields.keys():
            cls.n_fields[field] = dct_max_n_value[field]

    @property
    def mag_values(self):
        mag_values = OrderedDict([(filt, OrderedDict([(key, value) for key, value in self.items()
                                                      if re.findall(r'{}_mag_\d+$'.format(filt), key)]))
                                  for filt in StarFromCSV.filts])
        return mag_values

    def mag_digit_values(self, filt):
        return filter(lambda x: x != '-', tuple(self.mag_values[filt].values()))

    def reduction_AC(self, ABC_by_filt):
        filt = StarFromCSV.filts[0]
        a = ABC_by_filt[filt].a
        c = ABC_by_filt[filt].c

        for key, value in self.mag_values[filt].items():
            if value != '-':
                self[key] = a * value + c

        self.update_stat_mag_by_filter(filt)

    def reduction_ABC(self, ABC_by_filt):
        mag_values = copy.deepcopy(self.mag_values)
        for filt1 in StarFromCSV.filts:
            a = ABC_by_filt[filt1].a
            b = ABC_by_filt[filt1].b
            c = ABC_by_filt[filt1].c
            filt2 = self.choose_filter(filt1)
            for key_filt1, key_filt2 in zip(mag_values[filt1], mag_values[filt2]):
                if self[key_filt1] != '-' and self[key_filt1] != '-':
                    self[key_filt1] = a * mag_values[filt1][key_filt1] + b * (mag_values[filt2][key_filt2] - mag_values[filt1][key_filt1]) + c
                else:
                    self[key_filt1] = '-'

            self.update_stat_mag_by_filter(filt1)

    def mag_value_by_filter(self, filt):
        try:
            output = self.mag_values[filt]
        except KeyError as e:
            print("Filter's name error. You input {} filter's name ".format(e))
        else:
            return output

    def change_mag(self, change_lst, filt):
        for field, change in zip(self.mag_value_by_filter(filt).keys(), change_lst):
            if self[field] != '-':
                self[field] = float(self[field]) - change

        self.update_stat_mag_by_filter(filt)

    def update_stat_mag_by_filter(self, filt):
        mags = tuple(filter(lambda x: x != '-', self.mag_digit_values(filt)))
        if mags:
            self['{}_min_mag'.format(filt)] = max(mags)
            self['{}_max_mag'.format(filt)] = min(mags)
            self['{}_avr_mag'.format(filt)] = np.mean(mags)
            self['{}_std_mag'.format(filt)] = np.std(mags)
        else:
            self['{}_min_mag'.format(filt)] = '-'
            self['{}_max_mag'.format(filt)] = '-'
            self['{}_avr_mag'.format(filt)] = '-'
            self['{}_std_mag'.format(filt)] = '-'

    def update_mag(self):
        for filt in StarFromCSV.filts:
            self.update_stat_mag_by_filter(filt)

    @staticmethod
    def choose_filter(filt):
        if filt == 'V':
            filt2 = 'B'
        elif filt == 'R':
            filt2 = 'I'
        elif filt == 'B':
            filt2 = 'V'
        elif filt == 'I':
            filt2 = 'R'

        return filt2

# ----------------------------------------------------astrometry-------------------------------------------------------

def getheader(filepath):
    """
    Returns header of the file (filepath)
    :param filepath: path to file
    :return: header of the file (filepath)
    """
    with pyfits.open(filepath) as f:
        output = f[0].header

    return output

def writecsv(filepath, header, data):
    with open(filepath, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerows(data)

    print('File {} has created'.format(filepath))


def astrometry(dirpath):
    """
    Runs the solve-field soft for the list of fit-files. If header of fit-file has not empty values of 'RA' and 'DEC',
    then to the solve-field is transferred the values these fields. Renames new-files to fits-files (astrometry returns
     fits files with name: 'files.new'). Deletes excess files (astrometry returns many excess files).
    :param fit_path: path to fit-files
    :return: different files, that were produced during the operation, including fits-files with the extension new
    """
    fitfiles = [os.path.join(dirpath, filename) for filename in os.listdir(dirpath)
                if os.path.isfile(os.path.join(dirpath, filename))
                and filename.endswith('.fit')]
    solve_field = './solve-field'

    for filepath in sorted(fitfiles):
        header = getheader(filepath)
        try:
            if 'RA' in header and 'DEC' in header:
                if header['RA'] and header['DEC']:
                    os.system('{} --ra {} --dec {} --radius 0.5 --use-sextractor {}'.format(solve_field,
                                                                                        header['RA'],
                                                                                        header['DEC'],
                                                                                        filepath))
            else:
                os.system('{} --use-sextractor {}'.format(solve_field, filepath))
        except PermissionError as e:
            print(e)

        pre_filename = re.findall(r'^(.+)\.fit$', filepath)[0]
        try:
            os.rename('{}.{}'.format(pre_filename, 'new'),
                      '{}.{}'.format(pre_filename, 'fits'))
            os.system('chmod 777 {}'.format('{}.{}'.format(pre_filename, 'fits')))
        except OSError as e:
            print(e)

        try:
            os.remove(filepath)
        except OSError as e:
            print(e)

        delfiles = [os.path.join(dirpath, filename) for filename in os.listdir(dirpath)
                    if os.path.isfile(os.path.join(dirpath, filename))
                    and not filename.endswith('.fit')
                    and not filename.endswith('.fits')
                    and not filename.endswith('_cat.csv')]

        for filepath in delfiles:
            try:
                os.remove(filepath)
            except OSError as e:
                print(e)

# ---------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------S-Extractor-------------------------------------------------------


def sextractor(dirpath):
    """
    Runs s-extractor for files from the dirpath
    :param dirpath: the directory where ara placed fits-files
    :return: the result of the s-extractor
    """
    fitsfiles = [os.path.join(dirpath, filename) for filename in os.listdir(dirpath)
                 if filename.endswith('.fits')]

    m_date = OrderedDict()
    m_frame = list()

    for filepath in sorted(fitsfiles):

        file_frame, file_filt = re.findall(r'.*-(\d+)_(\w+)\.fits', filepath)[0]
        JD = str(getheader(filepath)['JD'])
        if file_frame not in m_frame:
            m_frame.append(file_frame)
        if file_filt in m_date:
            m_date[file_filt][file_frame] = JD
        else:
            m_date[file_filt] = {file_frame: JD}

        file_prename = re.findall(r'^(.+)\.fits$', filepath)[0]
        sexpath = '{}_sex.txt'.format(file_prename)
        try:
            os.system('sextractor -CATALOG_NAME {} -c daofind.sex {}'.format(sexpath, filepath))
        except OSError as e:
            print(e)
        except UnicodeDecodeError as e:
            print(e)

        data = []
        units = []
        pre_filename = re.findall(r'^(.+_sex)\.txt$', sexpath)[0]

        try:
            with open(sexpath, 'r') as f:
                for line in f:
                    if line[0] == '#':
                        units.append(line.strip().split()[2])
                    else:
                        data.append(line.strip().split())

            with open(sexpath, 'w', newline='') as f:
                f_writer = csv.writer(f, delimiter=',')
                f_writer.writerow(units)
                f_writer.writerows(data)


            os.rename('{}.{}'.format(pre_filename, 'txt'),
                      '{}.{}'.format(pre_filename, 'csv'))
            os.system('chmod 777 {}'.format('{}.{}'.format(pre_filename, 'csv')))
        except OSError as e:
            print(e)

    savepath = os.path.join(dirpath, 'frameJD.txt')
    with open(savepath, 'w') as f:
        f.write('Frame,{}\n'.format(','.join(m_date.keys())))
        for frame in sorted(m_frame):
            temp = ['{}'.format(frame)]
            for filt in m_date.keys():
                if frame in m_date[filt]:
                    temp.append(m_date[filt][frame])
                else:
                    temp.append('-')
            f.write('{}\n'.format(','.join(temp)))


# ---------------------------------------------------------------------------------------------------------------------

def inst_catalogue(dirpath, catpath):
    cross_with_catalogue(dirpath, catpath)
    fullinstcat(dirpath)


def column_from_csv(data, column_name):
    return list(map(lambda x: x[column_name], data))


def not_match_object(sexpath, XMatchfile):
    sex_header, sex_data = open_csv(sexpath, ',')
    match_header, match_data = open_csv(XMatchfile, '\t')
    coord = list(zip(column_from_csv(match_data, 'ALPHA_J2000_tab2'), column_from_csv(match_data, 'DELTA_J2000_tab2')))
    not_match_data = list(filter(lambda x: (x['ALPHA_J2000'], x['DELTA_J2000']) not in coord, sex_data))

    return not_match_data, sex_header


def update_catalogue(dirpath, sexpath, catpath):
    batpath = os.path.join(dirpath, 'bat.ajs')
    with open(batpath, 'w') as f:
        XMatchfile = '{}_XMatch.csv'.format(re.findall(r'(.+)_sex\.csv$', sexpath)[0])
        data_to_Aladin = ['#AJS',
                          'sync',
                          'load {}'.format(sexpath),
                          'load {}'.format(catpath),
                          'xmatch {} {} 1.0'.format(catpath.split(os.sep)[-1], sexpath.split(os.sep)[-1]),
                          'export XMatch {}'.format(XMatchfile),
                          'quit']

        for line in data_to_Aladin:
            f.write('{}\n'.format(line))

    os.system('java -jar /home/crao/Astronomy/Soft/Aladin/Aladin/Aladin.jar -script < {}'.format(batpath))
    # os.system('java -jar C:\Aladin.jar -script < {}'.format(batpath))

    try:
        not_match_data, sex_header = not_match_object(sexpath, XMatchfile)
        os.remove(XMatchfile)
    except OSError as e:
        print(e)

    cat_header, cat_data = open_csv(catpath, ',')
    for src in not_match_data:
        temp = {'RAJ2000': src['ALPHA_J2000'], 'DEJ2000': src['DELTA_J2000']}
        for src_1 in cat_header:
            if src_1 not in temp:
                temp[src_1] = '-'
        cat_data.append(temp)

    writecsv(catpath, cat_header, cat_data)




def cross_with_catalogue(dirpath, catpath):
    sexfiles = [os.path.join(dirpath, filename) for filename in os.listdir(dirpath)
                if filename.endswith('_sex.csv')]

    for sexpath in sorted(sexfiles):
        update_catalogue(dirpath, sexpath, catpath)
        batpath = os.path.join(dirpath, 'bat.ajs')
        with open(batpath, 'w') as f:
            XMatchfile = '{}_XMatch.csv'.format(re.findall(r'(.+)_sex\.csv$', sexpath)[0])
            data_to_Aladin = ['#AJS',
                              'sync',
                              'load {}'.format(sexpath),
                              'load {}'.format(catpath),
                              'xmatch {} {} 1.0'.format(catpath.split(os.sep)[-1], sexpath.split(os.sep)[-1]),
                              'export XMatch {}'.format(XMatchfile),
                              'quit']

            for line in data_to_Aladin:
                f.write('{}\n'.format(line))

        os.system('java -jar /home/crao/Astronomy/Soft/Aladin/Aladin/Aladin.jar -script < {}'.format(batpath))
        #os.system('java -jar C:\Aladin.jar -script < {}'.format(batpath))

        filt = re.findall(r'.*_(\w+)_sex.csv$', sexpath.split(os.sep)[-1])[0].upper()

        try:
            os.remove(batpath)

            with open(XMatchfile, 'r') as f:
                oldheader = f.readline().strip().split('\t')

            header = OrderedDict()
            for field in oldheader:
                if field == 'RAJ2000_tab1':
                    header[field] = 'RAJ2000'
                elif field == 'DEJ2000_tab1':
                    header[field] = 'DEJ2000'
                elif field == '{}_mag_tab1'.format(filt):
                    header[field] = '{}_mag_cat'.format(filt)
                elif re.findall(r'^MAG_.*_tab2$', field):
                    header[field] = '{}_mag_{}'.format(filt, re.findall(r'.*-(\d+)_\w+_sex.csv$', sexpath.split(os.sep)[-1])[0])

            with open(XMatchfile, 'r') as csvfile:
                csvfile = csv.DictReader(csvfile, delimiter='\t')
                next(csvfile)
                data = [dict(map(lambda key: (header[key], row[key]) if row[key] != ' ' else
                                                (header[key], '-'),
                                    filter(lambda x: x in header, list(row.keys()))))
                        for row in csvfile]

            writecsv(XMatchfile, list(header.values()), data)

            os.remove(sexpath)
        except OSError as e:
            print(e)

def getspace(temp_src, header):
    new_data = copy.deepcopy(temp_src)
    for field in header:
        if field not in new_data:
            new_data[field] = '-'
    return new_data


def getstat(data, m_data, filt0):
    for src in data:
        m = list(map(lambda key: float(src[key]), filter(lambda i: (i in m_data) and (src[i] != '-'), src.keys())))
        if m:
            src['{}_min_mag'.format(filt0)] = np.max(m)
            src['{}_max_mag'.format(filt0)] = np.min(m)
            src['{}_avr_mag'.format(filt0)] = np.mean(m)
            src['{}_std_mag'.format(filt0)] = np.std(m)
            src['{}_n'.format(filt0)] = len(m)
        else:
            src['{}_min_mag'.format(filt0)] = '-'
            src['{}_max_mag'.format(filt0)] = '-'
            src['{}_avr_mag'.format(filt0)] = '-'
            src['{}_std_mag'.format(filt0)] = '-'
            src['{}_n'.format(filt0)] = 0


def fullinstcat(dirpath):
    """
    Run cross-analysis for different filters
    :param dirpath: the directory where ara placed *_sex.txt-files
    :return: the files with data of all objects in all filters
    """

    matchfiles = [os.path.join(dirpath, filename) for filename in os.listdir(dirpath)
                 if filename.endswith('XMatch.csv')]

    if matchfiles:
        filt0 = re.findall(r'.*_(\w+)_XMatch.csv$',
                           sorted(matchfiles,
                                  key=lambda x: (re.findall(r'.*_(\w+)_XMatch.csv$', x)[0],
                                                 re.findall(r'.*-(\d+)_\w+_XMatch.csv$', x)[0]))[0])[0]
    header = ['RAJ2000', 'DEJ2000']
    m_data = []
    data = list()
    for index, filepath in enumerate(sorted(matchfiles, key=lambda x: (re.findall(r'.*_(\w+)_XMatch.csv$', x)[0],
                                                                       re.findall(r'.*-(\d+)_\w+_XMatch.csv$', x)[0]))):
        filt = re.findall(r'.*_(\w+)_XMatch.csv$', filepath)[0]

        if filt != filt0 or index == 0:
            getstat(data, m_data, filt0)
            m_data = []
            filt0 = filt

            if '{}_mag_cat'.format(filt0) not in header:
                header.append('{}_mag_cat'.format(filt0))
                header.append('{}_min_mag'.format(filt0))
                header.append('{}_max_mag'.format(filt0))
                header.append('{}_avr_mag'.format(filt0))
                header.append('{}_std_mag'.format(filt0))
                header.append('{}_n'.format(filt0))

            for src in data:
                if '{}_mag_cat'.format(filt0) not in src:
                    src['{}_mag_cat'.format(filt0)] = '-'

        with open(filepath, 'r') as csvfile:
            csvfile = csv.DictReader(csvfile, delimiter=',')
            temp_data = [row for row in csvfile]

        magfile = '{}_mag_{}'.format(filt0, re.findall(r'.*-(\d+)_\w+_XMatch.csv$', filepath)[0])

        header.append(magfile)
        m_data.append(magfile)

        if not data:
            data = [getspace(temp_src, header) for temp_src in temp_data]
        else:
            for src in data:
                src[magfile] = '-'

            for temp_src in temp_data:
                for src in data:
                    if src['RAJ2000'] == temp_src['RAJ2000'] and src['DEJ2000'] == temp_src['DEJ2000']:
                        src[magfile] = temp_src[magfile]
                        src['{}_mag_cat'.format(filt0)] = temp_src['{}_mag_cat'.format(filt0)]
                        break
                else:
                    data.append(getspace(temp_src, header))
        if index == len(matchfiles) - 1:
            getstat(data, m_data, filt0)
            n_field = list(filter(lambda x: re.findall(r'^\w+_n$', x), header))
            new_data = []
            for row in data:
                for field in n_field:
                    if row[field] == '-':
                        row[field] = 0
                for field in n_field:
                    if int(row[field]) > 1:
                        new_data.append(row)
                        break
        try:
            os.remove(filepath)
        except OSError as e:
            print(e)

    print('data = {}\nnew_data = {}\n'.format(len(data), len(new_data)))
    writecsv(os.path.join(dirpath,'{}_fullinstcat.csv'.format(re.findall(r'^(\w+?)_.*', matchfiles[0].split(os.sep)[-1])[0])),
             header, new_data)



def open_csv(csv_file_path, delimiter=','):
    with open(csv_file_path, 'r') as f:
        header = f.readline().strip().split(',')

    with open(csv_file_path, 'r') as csv_file:
        csv_file = csv.DictReader(csv_file, delimiter=delimiter)
        data = [row for row in csv_file]

    return header, data


def get_column(data, column_name):
    return list(map(lambda x: x[column_name], data))


def load_catalogue(dir_path):
    full_catalogue_path = [file.path for file in os.scandir(dir_path) if file.name.endswith('_fullinstcat.csv')][0]
    prefix_full_cat = re.findall(r'(.*)_fullinstcat.csv$', full_catalogue_path)[0]

    header, raw_data = open_csv(full_catalogue_path)
    data = convert_data_to_float(raw_data, header)

    return data, prefix_full_cat


def convert_data_to_float(raw_data, header):
    data = []
    unchange_field = ['RAJ2000', 'DEJ2000']
    for star in raw_data:
        temp = list(map(lambda field: (field, star[field]), unchange_field))
        temp.extend(list(map(lambda field: (field, float(star[field])) if star[field] != '-' else
                                           (field, '-'), filter(lambda field: field not in unchange_field, header))))
        data.append(StarFromCSV(temp))

    dct_max_n_value = dict([(field, max(data, key=lambda x: x[field])[field])
                            for field in StarFromCSV.n_fields.keys()])
    StarFromCSV.fill_n_fields(dct_max_n_value)

    return data


def build_catalogues(data, key):
    type_cat = dict(good=operator.add, to_sys=operator.mul)

    output_data_1 = []

    output_data = [src for src in data if reduce(type_cat[key],
                                               map(lambda field: src[field] == StarFromCSV.n_fields[field],
                                                   StarFromCSV.n_fields.keys()))]

    if key == 'to_sys':
        for src in output_data:
            temp = dict(map(lambda x: (x, src[x]),filter(lambda x: re.findall(r'.*_cat', x),  StarFromCSV.header))) #re.findall(r'*_mag_cat$', x)
            for key, value in temp.items():
                if value == '-':
                    break
            else:
                output_data_1.append(src)
    else:
        output_data_1 = output_data

    return copy.deepcopy(output_data_1)


def write_catalogue(dir_path, data, prefix, postfix, key):
    dct_file_name = dict(full='{}_full_{}.csv'.format(prefix, postfix),
                         good='{}_good_{}.csv'.format(prefix, postfix),
                         to_sys='{}_tosys_{}.csv'.format(prefix, postfix),
                         result='{}_result_{}.csv'.format(prefix, postfix))

    save_file_path = os.path.join(dir_path, dct_file_name[key])
    write_csv(save_file_path, data, StarFromCSV.header)


def write_csv(save_file_path, data, header):
    with open(save_file_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerows(data)

    if os.path.exists(save_file_path):
        print('File {} has created'.format(save_file_path))
    else:
        print('File {} has not created'.format(save_file_path))


def produce_corr_analysis(dir_path, data, to_sys_data, obj):
    ref_stars = ''

    #to_sys_data = smoothing(to_sys_data, 3)

    for filt in StarFromCSV.filts:
        temp_ten = sorted(filter(lambda x: '{}_{}'.format(x['RAJ2000'], x['DEJ2000']) != obj, to_sys_data),
                          key=lambda x: x['{}_std_mag'.format(filt)])[:10]

        ten_good = dict()
        tensns = dict()

        for index, src in enumerate(temp_ten):
            ten_good['{}_{}'.format(src['RAJ2000'], src['DEJ2000'])] = list(src.mag_digit_values(filt))
            tensns[index] = list(src.mag_digit_values(filt))

        rho = DataFrame(ten_good).corr()
        rhosns = DataFrame(tensns).corr()
        save_correlogram(dir_path, rhosns, filt)

        mycorr = get_max_corr(rho)

        ref1_RADE, ref2_RADE, err_arr = get_refstar(dir_path, to_sys_data, filt, mycorr, ten_good)

        changedata(to_sys_data, filt, err_arr)
        changedata(data, filt, err_arr)

    return ref_stars, data, to_sys_data


def get_refstar(dir_path, to_sys_data, filt, mycorr, ten_good):
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

    for mkey1, mkey2, maxr in sorted(mycorr, key=lambda x: x[2], reverse=True)[:5]:
        avr = np.mean(ten_good[mkey1])
        err_arr_0 = list(map(lambda x: x - avr, ten_good[mkey1]))

        temp_data = copy.deepcopy(to_sys_data)

        changedata(temp_data, filt, err_arr_0)
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

    return ref1_RADE, ref2_RADE, err_arr


def save_correlogram(dir_path, rhosns, filt):
    heatmap(rhosns, annot=True, fmt=".3f", linewidths=.5)
    plt.yticks(rotation=0)
    plt.title('Correlogram of stars with the least error ({} filter)\n'.format(filt), fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    savename = os.path.join(dir_path, 'corr_{}'.format(filt))
    plt.savefig(savename,
                  dpi=None,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None,
                  fmt='svg')
    plt.close()


def changedata(data, filt, err_arr):
    for star in data:
        star.change_mag(err_arr, filt)


def get_max_corr(rho):
    temp = []
    for key1 in rho:
        for key2 in rho:
            if key1 != key2:
                if temp:
                    for src in temp:
                        if src[0] == key2 and src[1] == key1:
                            break
                        elif src[0] == key1 and src[1] == key2:
                            break
                    else:
                        temp.append((key1, key2, rho[key1][key2]))
                else:
                    temp.append((key1, key2, rho[key1][key2]))

    temp.sort(key=lambda x: x[2], reverse=True)
    return temp


def get_mean_std(data, filt):
    stat = [src['{}_std_mag'.format(filt)] for src in data]
    mean = np.mean(stat)
    std = np.std(stat)
    return mean, std


def smoothing(to_sys_data, n=3):
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


def get_ABC(to_sys_data_corr):
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
    filt1_inst = []
    f2_f1_inst = []
    filt1_cat = []

    for src in data:
        filt1_inst.append(src['{}_avr_mag'.format(filt1)])
        f2_f1_inst.append(src['{}_avr_mag'.format(filt2)] - src['{}_avr_mag'.format(filt1)])
        filt1_cat.append(src['{}_mag_cat'.format(filt1)])

    return filt1_inst, f2_f1_inst, filt1_cat


def save_regress_line(dir_path, reducted_catalogue, ABC_by_filt):
    for filt in StarFromCSV.filts:
        x = get_column(reducted_catalogue, '{}_avr_mag'.format(filt))
        y = get_column(reducted_catalogue, '{}_mag_cat'.format(filt))

        A = np.vstack([x, np.ones(len(x))]).T
        k, b = np.linalg.lstsq(A, y)[0]

        y_mnk = [k * i + b for i in x]

        save_name = os.path.join(dir_path, 'regress_{}'.format(filt))
        xlab = 'Instrumental catalogue, mag'
        ylab = 'Catalogue, mag'
        tlt = 'Catalohue vs instrumental magnitudes ({} filter)'.format(filt)
        if ABC_by_filt[filt].b != 0:
            filt2 = StarFromCSV.choose_filter(filt)
            legendtxt = '${}_c={:.3}*{}_i{:+.3}*({}_i-{}_i){:+.3}$\n$y={:.3}*x{:+.3}$'.format(filt, ABC_by_filt[filt].a, filt, ABC_by_filt[filt].b, filt2, filt, ABC_by_filt[filt].c, k, b)
        else:
            legendtxt = '${}_c={:.3}*{}_i{:+.3}$\n$y={:.3}*x{:+.3}$'.format(filt, ABC_by_filt[filt].a, filt, ABC_by_filt[filt].c, k, b)
        save_mnk(x, y, y_mnk, legendtxt, xlab, ylab, tlt, save_name)

def save_mnk(x, y, y_mnk, legendtxt, xlab, ylab, tlt, savename):
    plt.plot(x, y_mnk, label=legendtxt)
    plt.legend(loc=2)
    plt.plot(x, y, 'ro')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(tlt)
    plt.grid(True)
    plt.savefig(savename,
                  dpi=None,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None)
    plt.close()


def save_regress(save_path, x, y, y_mnk, tlt, x_lbl, y_lbl, legend_lbl):
    plt.plot(x, y, linestyle=' ', marker='.', color='grey')
    plt.plot(x, y_mnk, color='red', label=legend_lbl)

    legend1 = plt.legend(loc=2, fontsize=16, prop={'size':14}, frameon=True, framealpha=0.5)
    frame = legend1.get_frame()
    frame.set_color('white')
    frame.set_edgecolor('black')

    plt.axis([9, 16, 9, 16])
    formatter = matplotlib.ticker.FormatStrFormatter("%.0f")
    plt.xticks(fontsize=14, fontname='Times')
    plt.yticks(fontsize=14, fontname='Times')
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().set_ylabel(y_lbl, fontsize=14, fontname='Times')
    plt.gca().set_xlabel(x_lbl, fontsize=14, fontname='Times')
    plt.gca().set_title(tlt, fontsize=16, fontname='Times')
    plt.grid(color='grey', linestyle='--', linewidth=0.5)

    plt.savefig(save_path,
                  dpi=None,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None,
                  fmt='svg')
    plt.close()


def make_reducted_catalogue(ABC_by_filt, data):
    make_data = copy.deepcopy(data)

    if len(StarFromCSV.filts) == 1:

        for star in make_data:
            star.reduction_AC(ABC_by_filt)


    else:
        for star in make_data:
            star.reduction_ABC(ABC_by_filt)

    return make_data


def load_JSON_list(path_to_json):
    with open(path_to_json) as json_file:
        inf_json = json.load(json_file)

    return inf_json


def go(dir_path):
    obj = '303.381794_+65.162131'

    #astrometry(dir_path)
    sextractor(dir_path)
    cat_path = [os.path.join(dir_path, file_name) for file_name in os.listdir(dir_path) if file_name.endswith('_cat.csv')][0]
    inst_catalogue(dir_path, cat_path)

    data, prefix_full_cat = load_catalogue(dir_path)
    to_sys_data = build_catalogues(data, 'to_sys')
    write_catalogue(dir_path, to_sys_data, prefix_full_cat, 'inst', 'to_sys')

    ref_stars, data_corr, to_sys_data_corr = produce_corr_analysis(dir_path, data, to_sys_data, obj)
    write_catalogue(dir_path, data_corr, prefix_full_cat, 'corr', 'full')
    write_catalogue(dir_path, to_sys_data_corr, prefix_full_cat, 'corr', 'to_sys')

    ABC_by_filt = get_ABC(to_sys_data_corr)
    reducted_catalogue = make_reducted_catalogue(ABC_by_filt, to_sys_data_corr)

    write_catalogue(dir_path, reducted_catalogue, prefix_full_cat, '', 'result')
    save_regress_line(dir_path, reducted_catalogue, ABC_by_filt)



def main():
    path_to_data = os.path.abspath(input('Enter the path to json:\n'))
    lst_json = load_JSON_list(path_to_data)
    #dir_path = os.path.abspath(input('Enter the path to directory:\n'))
    #obj = input('Enter the RA and DE of object (RA_DE):\n')
    for dir_path in lst_json:
        go(dir_path)
        StarFromCSV.reset_class_data()


def main1():
    dir_path = os.path.abspath(input('Enter the path to directory:\n'))
    go(dir_path)


if __name__ == '__main__':
    main1()
