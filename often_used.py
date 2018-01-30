import csv


def open_csv(csv_file_path, delimiter=',', first_line=False):
    '''
    Open csv file
    :param csv_file_path: Path to csv-file
    :param delimiter: delimiter in csv-file
    :param first_line: Skip (True) the first line or not (False)
    :return: header and data from csv-file
    '''

    with open(csv_file_path, 'r') as f:
        header = f.readline().strip().split(delimiter)

    with open(csv_file_path, 'r') as csv_file:
        csv_file = csv.DictReader(csv_file, delimiter=delimiter)
        if first_line:
            next(csv_file)
        stars_from_file = [row for row in csv_file]

    return header, stars_from_file


def write_to_csv(save_path, data, header, delimiter=','):
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
    print('File {} has created!'.format(save_path))


def replace_data_in_dct(dct, first_symbol, second_symbol):
    '''
    Replaces first_symbol with second_symbol in dictionary.
    :param dct: data from dictionary
    :param first_symbol: first symbol
    :param second_symbol: second symbol
    :return: data with replaced symbols
    '''

    for key, value in dct.items():
        if value == first_symbol:
            dct[key] = second_symbol

    return dct


def replace_data_in_dcts(lst, first_symbol, second_symbol):
    '''
    Replaces first_symbol with second_symbol in list with dictionaries.
    :param dct: list with dictionaries
    :param first_symbol: first symbol
    :param second_symbol: second symbol
    :return: data with replaced symbols
    '''

    lst = [replace_data_in_dct(dct, first_symbol, second_symbol) for dct in lst]

    return lst