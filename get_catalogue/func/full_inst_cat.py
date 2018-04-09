import os
import re
import numpy as np
from collections import OrderedDict
from func.work_with_csv import open_dict_csv, write_dict_csv, get_column


MAG_FIELD = '{}_mag_{}'


def full_inst_cat(dir_path, JD_path):
    '''
    Reades the processed cross-match catalogs between nomus and SExtractor. Takes the information about objects in this catalogs. 
    Writes the information in compilation catalogue. Saves the catalogue by the name name_fullinscat.csv.
    :param dir_path: path to folder with cross-match catalogs between nomus and SExtractor
    :param JD_path: path to file with information about the Julian Date in frames divided by filters.
    :return: path to compilation catalogue with information about all objects in cross-match catalogs between nomus and SExtractor
    '''
    fields, dct_fields, filts = fill_header_full_cat(JD_path)
    full_data, prefix_name_cat = fill_full_cat(dir_path, fields)
    calc_stat_data(full_data, dct_fields, filts)
    save_path = os.path.join(dir_path, '{}_fullinstcat.csv'.format(prefix_name_cat))
    write_dict_csv(save_path, fields, full_data)
    return save_path


def fill_header_full_cat(JD_path):
    '''
    Reades the file with file with information about the Julian Date in frames divided by filters. 
    Convert the information to header of the compilation catalogue.
    :param JD_path: path to file with information about the Julian Date in frames divided by filters.
    :return: header and header divided by filters of the compilation catalogue, list of filters placed in the compilation catalogue.
    '''
    JD_header, JD_data = open_dict_csv(JD_path)
    filts = list(filter(lambda x: x != 'Frame', JD_header))

    stat_fields = ['{}_mag_cat', '{}_min_mag', '{}_max_mag', '{}_avr_mag', '{}_std_mag', '{}_n']
    dct_fields = OrderedDict([(filt, [src.format(filt) for src in stat_fields]) for filt in filts])
    for src in JD_data:
        for filt in filts:
            if src[filt] != '-':
                dct_fields[filt].append(MAG_FIELD.format(filt, src['Frame']))
    fields = ['RAJ2000', 'DEJ2000']
    for value in dct_fields.values():
        fields.extend(value)

    return fields, dct_fields, filts 


def fill_full_cat(dir_path, fields):
    '''
    Reades the processed cross-match catalogs between nomus and SExtractor. Takes the information about objects in this catalogs. 
    Writes the information in the compilation catalogue.
    :param dir_path: path to folder with cross-match catalogs between nomus and SExtractor.
    :param fields: header of the compilation catalogue.
    :return: the compilation catalogue without statistical information about objects and prefix of name of catalogue.
    '''
    CAT_MAG_FIELD = '{}_mag_cat'
    NUMBER_FILTER_XM = r'.*-(\d+)_(\w+)_XMatch.csv$'
    match_files = [file for file in os.scandir(dir_path) if file.name.endswith('XMatch.csv')]
    match_files.sort(key=lambda x: re.findall(NUMBER_FILTER_XM, x.name)[0][::-1])
    full_data = list()
    for file_XM in match_files:
        m_num, filt = re.findall(NUMBER_FILTER_XM, file_XM.name)[0]
        XM_header, XM_data = open_dict_csv(file_XM.path)
        XM_mag_field = MAG_FIELD.format(filt, m_num)
        XM_cat_mag_field = CAT_MAG_FIELD.format(filt)
        for star in XM_data:
            if full_data:
                stars_in_full_cat = list(zip(get_column(full_data, 'RAJ2000'), get_column(full_data, 'DEJ2000')))
                star_coord = (star['RAJ2000'], star['DEJ2000'])
                if star_coord in stars_in_full_cat:
                    temp_index = stars_in_full_cat.index(star_coord)
                    full_data[temp_index][XM_mag_field] = str(round(float(star[XM_mag_field]), 4))
                    full_data[temp_index][XM_cat_mag_field] = star[XM_cat_mag_field]
                else:
                    full_data.append(OrderedDict([(field, '-') for field in fields]))
                    full_data[-1]['RAJ2000'] = star['RAJ2000']
                    full_data[-1]['DEJ2000'] = star['DEJ2000']
                    full_data[-1][XM_mag_field] = str(round(float(star[XM_mag_field]), 4))
                    full_data[-1][XM_cat_mag_field] = star[XM_cat_mag_field]
            else:
                full_data.append(OrderedDict([(field, '-') for field in fields]))
                full_data[-1]['RAJ2000'] = star['RAJ2000']
                full_data[-1]['DEJ2000'] = star['DEJ2000']
                full_data[-1][XM_mag_field] = str(round(float(star[XM_mag_field]), 4))
                full_data[-1][XM_cat_mag_field] = star[XM_cat_mag_field]
    full_data.sort(key=lambda x: x['RAJ2000'])
    prefix_name_cat = re.findall(r'^(\w+?)-.*', match_files[0].name)[0]

    return full_data, prefix_name_cat



def calc_stat_data(full_data, dct_fields, filts):
    '''
    Computes statistical information about objects. Adds this information to the compilation catalogue.
    :param full_data: information about object in the compilation catalogue.
    :param dct_fields: header divided by filters of the compilation catalogue.
    :param filts: list of filters placed in the compilation catalogue.
    :return: statistical information about objects. Writes it to the compilation catalogue.
    '''
    for star in full_data:
        for filt in filts:
            m_fields = list(filter(lambda x: star[x] != '-', filter(lambda x: re.findall(r'{}_mag_\d+'.format(filt), x), dct_fields[filt])))
            value_m_fields = list(map(lambda x: float(star[x]),m_fields))
            if value_m_fields:
                star['{}_min_mag'.format(filt)] = min(value_m_fields)
                star['{}_max_mag'.format(filt)] = max(value_m_fields)
                star['{}_avr_mag'.format(filt)] = round(np.mean(value_m_fields), 4)
                star['{}_std_mag'.format(filt)] = round(np.std(value_m_fields), 4)
                star['{}_n'.format(filt)] = len(value_m_fields)
            else:
                star['{}_min_mag'.format(filt)] = '-'
                star['{}_max_mag'.format(filt)] = '-'
                star['{}_avr_mag'.format(filt)] = '-'
                star['{}_std_mag'.format(filt)] = '-'
                star['{}_n'.format(filt)] = '0'
