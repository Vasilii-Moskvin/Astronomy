from collections import OrderedDict, namedtuple
import csv
import os.path
import math


ALADIN_PATH = r'C:\Aladin.jar'


def open_csv(csv_file_path, delimiter=',', first_line=False):
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
    if os.path.exists(save_path):
        print('File {} has created!'.format(save_path))
    else:
        print('Error: File {} has no created!'.format(save_path))


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

def column_from_list_dct(data, dct_key):
    '''
    Gets a list of data from the list of dictionaries by key in dictionaries.
    :param data: list of dictionaries
    :param dct_key: key in dictionaries
    :return: list of data on the key in dictionaries.
    '''
    return list(map(lambda x: x[dct_key], data))


def get_XMatch_NOMAD_USNO(save_dir_path, ra, de, radius="14'"):
    """
    Gets XMatch-file from Aladin.
    :param save_dir_path: path to save directory for catalogue from Aladin
    :param ra: RA
    :param de: DE
    :param radius: radius of area for catalogue
    :return: path to XMatch-file from Aladin
    """
    global ALADIN_PATH
    XMatch_NOMAD_USNO_path = os.path.join(save_dir_path, 'nomus_cat.csv')
    NOMAD_path = os.path.join(save_dir_path, 'NOMAD_nomus_cat.csv')
    USNO_path = os.path.join(save_dir_path, 'USNO_nomus_cat.csv')
    script_for_Aladin = ['#AJS',
                         'sync',
                         'get VizieR(NOMAD1) {} {} {}'.format(ra, de, radius),
                         'get VizieR(USNO-B1) {} {} {}'.format(ra, de, radius),
                         'xmatch NOMAD1 USNO-B1 1.0',
                         'export NOMAD1 {}'.format(NOMAD_path),
                         'export USNO-B1 {}'.format(USNO_path),
                         'export XMatch {}'.format(XMatch_NOMAD_USNO_path),
                         'quit']

    script_for_Aladin_path = os.path.join(save_dir_path, 'script_for_Aladin.txt')
    with open(script_for_Aladin_path, 'w') as w:
        for line in script_for_Aladin:
            w.write('{}\n'.format(line))

    os.system('java -jar {} -script < {}'.format(ALADIN_PATH, script_for_Aladin_path))
    os.remove(script_for_Aladin_path)

    if os.path.exists(XMatch_NOMAD_USNO_path):
        return XMatch_NOMAD_USNO_path, NOMAD_path, USNO_path
    else:
        return False


def convert_XMatch_to_catalogue( stars_inform, epoch=''):
    """
    Converts XMatch-file from Aladin to csv- catalogue.
    :param NOMAD_USNO_path: path to XMatch NOMAD1 and USNO-B1 catalogues
    :return: csv- catalogue
    """
    if epoch:
        convert_coord_NOMAD = OrderedDict((('_RAJ{}_tab1'.format(epoch), 'RAJ{}'.format(epoch)),
                                           ('_DEJ{}_tab1'.format(epoch), 'DEJ{}'.format(epoch))))
        convert_coord_USNO = OrderedDict((('RAJ{}_tab2'.format(epoch), 'RAJ{}'.format(epoch)),
                                          ('DEJ{}_tab2'.format(epoch), 'DEJ{}'.format(epoch))))
        convert_data = OrderedDict((('NOMAD1_tab1', 'NOMAD1_ID'),
                                      ('USNO-B1.0_tab2', 'USNO-B1.0'),
                                      ('Bmag_tab1', 'B_mag'),
                                      ('Vmag_tab1', 'V_mag'),
                                      ('Rmag_tab1', 'R_mag'),
                                      ('Imag_tab2', 'I_mag')))
    else:
        convert_coord_NOMAD = OrderedDict((('_RAJ2000_tab1', 'RAJ2000'),
                                           ('_DEJ2000_tab1', 'DEJ2000')))
        convert_coord_USNO = OrderedDict((('RAJ2000_tab2', 'RAJ2000'),
                                          ('DEJ2000_tab2', 'DEJ2000')))
        convert_data = OrderedDict((('NOMAD1_tab1', 'NOMAD1_ID'),
                                      ('USNO-B1.0_tab2', 'USNO-B1.0'),
                                      ('Bmag_tab1', 'B_mag'),
                                      ('Vmag_tab1', 'V_mag'),
                                      ('Rmag_tab1', 'R_mag'),
                                      ('Imag_tab2', 'I_mag')))
    save_header = list(convert_coord_NOMAD.values())
    save_header.extend(list(convert_data.values()))
    data_to_csv = []
    for star in stars_inform:
        convert_header = OrderedDict()
        if star['_RAJ2000_tab1'] != '-':
            convert_header.update(convert_coord_NOMAD)
        else:
            convert_header.update(convert_coord_USNO)
        convert_header.update(convert_data)
        temp_dct = dict(map(lambda key: (convert_header[key], star[key]),
                            filter(lambda x: x in convert_header, star.keys())))
        data_to_csv.append(temp_dct)

    return data_to_csv, save_header


def calc_today_coordinates(stars_from_file, epoch):
    '''
    Calculates the objects' coordinates, which are in the XMatch_NOMAD_USNO_path file for the epoch,
    indicated in the variable epoch.
    :param stars_from_file: XMatch NOMAD1 and USNO-B1 catalogues
    :param epoch: epoch
    :return: XMatch NOMAD1 and USNO-B1 catalogues with added colums of new coordinates
    '''

    float_epoch = float(epoch) - 2000

    for star in stars_from_file:
        if star['_RAJ2000_tab1'] != '-':
            ra = '_RAJ2000_tab1'
            dec = '_DEJ2000_tab1'
            pm_ra = 'pmRA_tab1'
            pm_dec = 'pmDE_tab1'
            ra_by_epoch = '_RAJ{}_tab1'.format(epoch)
            dec_by_epoch = '_DEJ{}_tab1'.format(epoch)
        else:
            ra = 'RAJ2000_tab2'
            dec = 'DEJ2000_tab2'
            pm_ra = 'pmRA_tab2'
            pm_dec = 'pmDE_tab2'
            ra_by_epoch = 'RAJ{}_tab2'.format(epoch)
            dec_by_epoch = 'DEJ{}_tab2'.format(epoch)


        star[dec_by_epoch] = round(float(star[dec]) + (0.001 * float(star[pm_dec]) * float_epoch) / 3600, 7)
        old_cos_dec = math.cos((math.pi / 180) * float(star[dec]))
        new_cos_dec = math.cos((math.pi / 180) * float(star[dec_by_epoch]))
        star[ra_by_epoch] = round(float(star[ra]) + \
                            ((0.001 * float(star[pm_ra]) * float_epoch) / 3600) * (new_cos_dec / old_cos_dec), 7)

    return stars_from_file


def add_all_stars(stars_XMatch, header_NOMAD1, stars_NOMAD1, header_USNO, stars_USNO):
    '''
    Adds stars to the XMatch catalogue for which no intersections have been found in the NOMAD1 and USNO-B1 catalogs.
    :param stars_XMatch: data about stars from XMatch catalogue
    :param header_NOMAD1: header of data from NOMAD1 catalogue
    :param stars_NOMAD1: data about stars from NOMAD1 catalogue
    :param header_USNO: header of data from USNO-B1 catalogue
    :param stars_USNO: data about stars from USNO-B1 catalogue
    :return:
    '''
    NOMAD_ID_in_XMatch = set(column_from_list_dct(stars_XMatch, 'NOMAD1_tab1'))
    USNO_ID_in_XMatch = set(column_from_list_dct(stars_XMatch, 'USNO-B1.0_tab2'))

    NOMAD_ID_in_NOMAD = set(column_from_list_dct(stars_NOMAD1, 'NOMAD1'))
    NOMAD_ID_in_USNO = set(column_from_list_dct(stars_USNO, 'USNO-B1.0'))

    delta_NOMAD = list(NOMAD_ID_in_NOMAD - NOMAD_ID_in_XMatch)
    delta_USNO = list(NOMAD_ID_in_USNO - USNO_ID_in_XMatch)

    NOMAD_not_in_XMatch = list(filter(lambda x: x['NOMAD1'] in delta_NOMAD, stars_NOMAD1))
    USNO_not_in_XMatch = list(filter(lambda x: x['USNO-B1.0'] in delta_USNO, stars_USNO))

    empty_NOMAD_in_XMatch = OrderedDict(map(lambda x: ('{}_tab1'.format(x), ' '), header_NOMAD1))
    empty_USNO_in_XMatch = OrderedDict(map(lambda x: ('{}_tab2'.format(x), ' '), header_USNO))

    convert_NOMAD_in_XMatch = OrderedDict(map(lambda x: (x, '{}_tab1'.format(x)), header_NOMAD1))
    convert_USNO_in_XMatch = OrderedDict(map(lambda x: (x, '{}_tab2'.format(x)), header_USNO))

    temp_lst = []
    for src in NOMAD_not_in_XMatch:
        temp_dct = OrderedDict(map(lambda x: (convert_NOMAD_in_XMatch[x], src[x]), src.keys()))
        temp_dct.update(empty_USNO_in_XMatch)
        temp_lst.append(temp_dct)

    for src in USNO_not_in_XMatch:
        temp_dct = OrderedDict()
        temp_dct.update(empty_NOMAD_in_XMatch)
        temp_dct.update(OrderedDict(map(lambda x: (convert_USNO_in_XMatch[x], src[x]), src.keys())))
        temp_lst.append(temp_dct)

    stars_XMatch.extend(temp_lst)

    return stars_XMatch



def main():
    save_dir_path = os.path.abspath(input('Enter the path where to save NomUs catalogue:\n'))
    epoch = input('Enter Epoch:\n')
    ra = input('Enter RA:\n')
    de = input('Enter DE:\n')
    radius = input('Enter radius:\n')

    try:
        XMatch_NOMAD_USNO_path, NOMAD_path, USNO_path = get_XMatch_NOMAD_USNO(save_dir_path, ra, de, radius)
        if XMatch_NOMAD_USNO_path and NOMAD_path and USNO_path:
            header_XMatch, stars_XMatch = open_csv(XMatch_NOMAD_USNO_path, delimiter='\t', first_line=True)
            header_NOMAD1, stars_NOMAD1 = open_csv(NOMAD_path, delimiter='\t', first_line=True)
            header_USNO, stars_USNO = open_csv(USNO_path, delimiter='\t', first_line=True)
            stars_all = add_all_stars(stars_XMatch, header_NOMAD1, stars_NOMAD1, header_USNO, stars_USNO)
            stars_inform = [replace_data_in_dct(dct, ' ', '-') for dct in stars_all]
            if epoch:
                stars_inform = calc_today_coordinates(stars_inform, epoch)
            stars_inform, save_header = convert_XMatch_to_catalogue(stars_inform, epoch)
            write_to_csv(XMatch_NOMAD_USNO_path, stars_inform, save_header, delimiter=',')
            os.remove(NOMAD_path)
            os.remove(USNO_path)
    except OSError as e:
        print(e)


if __name__ == '__main__':
    main()
