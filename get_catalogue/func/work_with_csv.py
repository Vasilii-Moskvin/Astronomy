import csv
import os.path

def open_dict_csv(csv_file_path, delimiter=',', first_line=False):
    '''
    Open csv file
    :param csv_file_path: Path to csv-file
    :param delimiter: delimiter in csv-file
    :param first_line: Skip (True) the first line or not (False)
    :return: header and data from csv-file
    '''
    with open(csv_file_path, 'r') as csv_file:
        csv_file = csv.DictReader(csv_file, delimiter=delimiter)
        header = csv_file.fieldnames
        if first_line:
            next(csv_file)
        data = [row for row in csv_file]

    return header, data


def write_dict_csv(save_path, header, data, delimiter=','):
    '''
    Writes data in csv-file.
    :param save_path: savefile path
    :param data: data
    :param header: header of data
    :param delimiter: delimiter
    :return: csv-file with data
    '''
    with open(save_path, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=header, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(data)
    if os.path.exists(save_path):
        print('File {} has created!'.format(save_path))
    else:
        print('Error: File {} has no created!'.format(save_path))


def write_csv(save_path, header, data, delimiter=','):
    '''
    Writes data in csv-file.
    :param save_path: savefile path
    :param data: data
    :param header: header of data
    :param delimiter: delimiter
    :return: csv-file with data
    '''
    with open(save_path, 'w', newline='') as csv_file:
        f_writer = csv.writer(csv_file, delimiter=delimiter)
        f_writer.writerow(header)
        f_writer.writerows(data)
    if os.path.exists(save_path):
        print('File {} has created!'.format(save_path))
    else:
        print('Error: File {} has no created!'.format(save_path))


def get_column(data, field_name):
    '''
    Gets a list of dictionaries and a field name from the dictionaries. Return a list of values from the dictionaries by the field name.
    :param data: a list of dictionaries.
    :param field_name: a field name from the dictionaries.
    :return: list of values from the dictionaries by the field name.
    '''
    return list(map(lambda x: x[field_name], data))
