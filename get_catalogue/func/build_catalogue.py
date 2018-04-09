import copy
import operator
import os
import re
from functools import reduce
from func.StarFromCSV import StarFromCSV
from func.work_with_csv import open_dict_csv, write_dict_csv


def load_catalogue(dir_path, full_inst_cat_path):
    '''
    Loads the compilation catalogue.
    :param dir_path: path to the folder with the compilation catalogue.
    :param full_inst_cat_path: path to the compilation catalogue.
    :return: Data from the compilation catalogue.
    '''
    prefix_full_cat = re.findall(r'(.*)_fullinstcat.csv$', full_inst_cat_path)[0]
    header, raw_data = open_dict_csv(full_inst_cat_path)
    data = convert_data_to_float(header, raw_data)

    return data, prefix_full_cat


def convert_data_to_float(header, raw_data):
    '''
    Converts data from the compilation catlogue to data consisting of floating-point numbers.
    :param header: header of the compilation catlogue.
    :param raw_data: raw data from the compilation catalogue.
    :return: data from the compilation catalogue consisting of floating-point numbers.
    '''
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
    '''
    Gets data that must be written to the catalogue.
    :param data: data from the compilation catalogue.
    :param key: type of recordable catalogue.
    :return: data that must be written to the catalogue.
    '''
    type_cat = dict(good=operator.add, to_sys=operator.mul)
    output_data_1 = []
    output_data = [src for src in data if reduce(type_cat[key],
                                                 map(lambda field: src[field] == StarFromCSV.n_fields[field],
                                                                   StarFromCSV.n_fields.keys()))]
    if key == 'to_sys':
        for src in output_data:
            temp = dict(map(lambda x: (x, src[x]),filter(lambda x: re.findall(r'.*_cat', x),  StarFromCSV.header)))
            for key, value in temp.items():
                if value == '-':
                    break
            else:
                output_data_1.append(src)
    else:
        output_data_1 = output_data

    return copy.deepcopy(output_data_1)


def write_catalogue(dir_path, data, prefix, postfix, key):
    '''
    Writes data to a catalogue. Constructs name of the catalogue by prefix and postfix. 
    :param dir_path: path to the folder with the compilation catalogue.
    :param data: data that must be written to the catalogue.
    :param prefix: the catalogue name prefix.
    :param postfix: the catalogue name postfix.
    :param key: type of recordable catalogue.
    :return: csv-files with the catalogue.
    '''
    dct_file_name = dict(full='{}_full_{}.csv'.format(prefix, postfix),
                         good='{}_good_{}.csv'.format(prefix, postfix),
                         to_sys='{}_tosys_{}.csv'.format(prefix, postfix),
                         result='{}_result_{}.csv'.format(prefix, postfix),)

    save_file_path = os.path.join(dir_path, dct_file_name[key])
    write_dict_csv(save_file_path, StarFromCSV.header, data)
