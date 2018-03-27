import os
import re
from func.work_with_csv import open_dict_csv, write_dict_csv, get_column
from collections import OrderedDict


PATH_TO_ALADIN = r'/home/crao/Astronomy/Soft/Aladin/Aladin/Aladin.jar'


def cross_with_catalogue(dir_path, cat_path, two_flag):
    '''
    Reads catalogs obtained from SExtractor. Using the Aladin software produces 
    a catalog of cross-matches between the nomus SExtractor and catalogs.
    :param dir_path: path to folder with SExtracor catalogs.
    :param cat_path: path to the reference catalogue.
    :param two_flag: if True then all objects in the image are searched,
                     and not only those that are in the reference directory.
    :return: saves cross-match catalogue in the folder dir_path.
    '''
    sex_files = [file.path for file in os.scandir(dir_path) if file.name.endswith('_sex.csv')]
    for sex_path in sorted(sex_files):
        if two_flag:
            update_catalogue(dirpath, sexpath, catpath)
        XMatchfile = get_XMatchcat_from_Aladin(dir_path, sex_path, cat_path)
        Aladin_cat_to_csv(sex_path, XMatchfile)


def get_XMatchcat_from_Aladin(dir_path, sex_path, cat_path):
    '''
    Using the Aladin software produces a catalog of cross-matches between the nomus SExtractor and catalogs.
    :param dir_path: path to folder with SExtracor catalogs.
    :param sex_path: path to SExtracor catalogue.
    :param cat_path: path to the reference catalogue.
    :return: name of the cross-match catalogue in the folder dir_path. 
    '''
    bat_path = os.path.join(dir_path, 'bat.ajs')
    XMatchfile = '{}_XMatch.csv'.format(re.findall(r'(.+)_sex\.csv$', sex_path)[0])
    data_to_Aladin = ['#AJS',
                      'sync',
                      'load {}'.format(sex_path),
                      'load {}'.format(cat_path),
                      'xmatch {} {} 1.0'.format(cat_path.split(os.sep)[-1], sex_path.split(os.sep)[-1]),
                      'export XMatch {}'.format(XMatchfile),
                      'quit']
    with open(bat_path, 'w') as f:
        for line in data_to_Aladin:
            f.write('{}\n'.format(line))

    os.system('java -jar {} -script < {}'.format(PATH_TO_ALADIN, bat_path))
    try:
        os.remove(bat_path)
    except OSError as e:
        print(e)
    
    return XMatchfile


def Aladin_cat_to_csv(sex_path, XMatchfile):
    '''
    Converts the cross-match catalogue obtained from Aladin software to csv format. Removes unwanted columns.
    Deletes the catalogue obtained from SExtracor.
    :param sex_path: path to catalogue obtained from SExtracor.
    :param cat_path: path to the reference catalogu.
    :return: saves cross-match catalogue in the folder dir_path.
    '''
    filt = re.findall(r'.*_(\w+)_sex.csv$', sex_path.split(os.sep)[-1])[0].upper()

    try:
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
                header[field] = '{}_mag_{}'.format(filt, re.findall(r'.*-(\d+)_\w+_sex.csv$', sex_path.split(os.sep)[-1])[0])
          
        XM_header, XM_data = open_dict_csv(XMatchfile, delimiter='\t', first_line=True)
        data = [dict(map(lambda key: (header[key], row[key]) if row[key] != ' ' else
                                     (header[key], '-'),
                                     filter(lambda x: x in header, list(row.keys())))) for row in XM_data]
        write_dict_csv(XMatchfile, list(header.values()), data)
        os.remove(sex_path)
    except OSError as e:
        print(e)


def not_match_object(sex_path, XMatchfile):
    '''
    Data about objects that are not in the reference catalogue.
    :param sex_path: path to catalogue obtained from SExtracor.
    :param XMatchfile: path to the cross-match catalogue between catalogue obtained from SExtracor and the reference catalogue.
    :return: data about objects that are not in the reference catalogue.
    '''
    sex_header, sex_data = open_dict_csv(sex_path, ',')
    match_header, match_data = open_dict_csv(XMatchfile, '\t')
    coord = list(zip(get_column(match_data, 'ALPHA_J2000_tab2'), get_column(match_data, 'DELTA_J2000_tab2')))
    not_match_data = list(filter(lambda x: (x['ALPHA_J2000'], x['DELTA_J2000']) not in coord, sex_data))

    return not_match_data, sex_header


def update_catalogue(dir_path, sex_path, cat_path):
    '''
    Adds to the reference catalogue the objects that are on the frame, but not in the reference catalogue.
    :param dir_path: path to work folder.
    :param sex_path: path to catalogue obtained from SExtracor.
    :param cat_path: path to the reference catalogu.
    :return: updated the reference catalogue. Adds to it the objects that are on the frame, but not in the reference catalogue.
    '''
    bat_path = os.path.join(dir_path, 'bat.ajs')
    with open(bat_path, 'w') as f:
        XMatchfile = '{}_XMatch.csv'.format(re.findall(r'(.+)_sex\.csv$', sex_path)[0])
        data_to_Aladin = ['#AJS',
                          'sync',
                          'load {}'.format(sex_path),
                          'load {}'.format(cat_path),
                          'xmatch {} {} 1.0'.format(cat_path.split(os.sep)[-1], sex_path.split(os.sep)[-1]),
                          'export XMatch {}'.format(XMatchfile),
                          'quit']

        for line in data_to_Aladin:
            f.write('{}\n'.format(line))

    os.system('java -jar /home/crao/Astronomy/Soft/Aladin/Aladin/Aladin.jar -script < {}'.format(bat_path))
    # os.system('java -jar C:\Aladin.jar -script < {}'.format(bat_path))

    try:
        not_match_data, sex_header = not_match_object(sex_path, XMatchfile)
        os.remove(XMatchfile)
    except OSError as e:
        print(e)

    cat_header, cat_data = open_dict_csv(cat_path, ',')
    for src in not_match_data:
        temp = {'RAJ2000': src['ALPHA_J2000'], 'DEJ2000': src['DELTA_J2000']}
        for src_1 in cat_header:
            if src_1 not in temp:
                temp[src_1] = '-'
        cat_data.append(temp)

    write_dict_csv(catpath, cat_header, cat_data)
