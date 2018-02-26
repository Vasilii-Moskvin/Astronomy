import os
import csv
from collections import OrderedDict
import numpy as np
import copy
import json


keep_keys = set(['JD', 'JD_0', 'GJD', 'GJD_0', 'phase', 'G_phase'])


def load_JD0_period(path_to_json):
    with open(path_to_json) as json_file:
        inf_json = json.load(json_file)

    null_epoch = inf_json['null_epoch']
    period = inf_json['period']

    return null_epoch, period


def column_from_csv(data, column_name):
    return list(map(lambda x: x[column_name], data))


def open_csv(file_path):
    with open(file_path, 'r') as f:
        header = f.readline().strip().split(',')

    with open(file_path, 'r') as csv_file:
        csv_file = csv.DictReader(csv_file, delimiter=',')
        raw_data = [row for row in csv_file]

    return header, raw_data


def get_mean_err(lst, n):
    mean_x = []
    err_x = []
    temp_x = []
    count = 1

    for x in lst:
        if count == n:
            count = 1
            mean_x.append(round(np.mean(temp_x), 6))
            err_x.append(round(np.std(temp_x), 6))
            temp_x = []
        temp_x.append(x)
        count += 1

    return mean_x, err_x



def add_error_to_data(data, n=3):
    copy_data = copy.deepcopy(data)
    for key, value in copy_data.items():
        new_value, err_value = get_mean_err(value, n)
        if key in keep_keys:
            data['{}_new'.format(key)] = new_value
        else:
            data['{}_new'.format(key)] = new_value
            data['{}_err'.format(key)] = err_value


def write_new_data(save_path, data, header):
    with open(save_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerows(data)
    print('File {} has created!'.format(save_path))


def transform_data_to_DictWriter(data):
    new_data = []
    for index, src in enumerate(data['JD']):
        new_data.append(OrderedDict(map(lambda field: (field, data[field][index])
                                        if len(data[field]) > index else (field, ''), list(data.keys()))))

    return new_data


def add_phase(data, header, JD_0, period):
    #JD_0 = 2457796.22269676
    #JD_0 = float(data[0]['JD'])
    for src in data:
        src['phase'] = ((float(src['JD']) - JD_0) % period) / period

    header.insert(2, 'phase')
    data.sort(key=lambda x: x['phase'])


def add_JD_0(data, header):
    JD_0 = int(data[0]['JD'].split('.')[0])
    for src in data:
        src['JD_0'] = (float(src['JD']) - JD_0)

    header.insert(1, 'JD_0')


def add_d_mag(data, header):
    n = 30
    d_keys = tuple(set(header) - keep_keys)

    dct_data = dict(map(lambda x: (x, column_from_csv(data, x)), d_keys))
    dct_mean_data = dict()
    for key, value in dct_data.items():
        temp = value[:n]
        temp.extend(value[:-n - 1:-1])
        dct_mean_data[key] = float(np.mean(list(map(float, temp))))

    for src in data:
        temp = copy.deepcopy(src)
        for key, value in temp.items():
            if key in d_keys:
                new_key = 'd{}'.format(key)
                src[new_key] = float(value) - dct_mean_data[key]

    for key in d_keys:
        new_key = 'd{}'.format(key)
        header.append(new_key)

def main_start(file_path, null_epoch, period):
    dir_path = os.sep.join(file_path.split(os.sep)[:-1])
    name_for_new_file = os.path.join(dir_path, '{}_new.csv'.format(file_path.split(os.sep)[-1].split('.')[0]))

    header, raw_data = open_csv(file_path)
    add_JD_0(raw_data, header)
    add_phase(raw_data, header, float(null_epoch), float(period))
    add_d_mag(raw_data, header)
    data = OrderedDict(map(lambda field: (field,
                                          list(map(lambda value_from_data: float(value_from_data[field]),
                                                   raw_data))), header))
    add_error_to_data(data, 10)
    new_data = transform_data_to_DictWriter(data)
    write_new_data(name_for_new_file, new_data, list(data.keys()))


def main():
    file_path = os.path.abspath(input('Enter the path to csv-file:\n'))
    path_to_json = os.path.abspath(input('Enter the path to json-file with JD0 and period:\n'))
    #path_to_json = 'C:\\Users\\vasil\\Desktop\\NSVS01031772\\Results\\DAT_NSVS01031772.json'
    null_epoch, period = load_JD0_period(path_to_json)

    main_start(file_path, null_epoch, period)


def main1():
    dir_path = os.path.abspath(input('Enter the path to directory with csv-files:\n'))
    for file in os.scandir(dir_path):
        file_path = file.path
        #file_path = os.path.abspath(input('Enter the path to csv-file:\n'))
        #path_to_json = os.path.abspath(input('Enter the path to json-file with JD0 and period:\n'))
        path_to_json = 'C:\\Users\\vasil\\Desktop\\NSVS01031772\\Results\\DAT_NSVS01031772.json'
        null_epoch, period = load_JD0_period(path_to_json)

        main_start(file_path, null_epoch, period)



if __name__ == '__main__':
    main()