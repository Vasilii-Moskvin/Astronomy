import os
import re
from func.work_with_csv import open_dict_csv, write_dict_csv
from collections import OrderedDict


PATH_TO_ALADIN = r'/home/crao/Astronomy/Soft/Aladin/Aladin/Aladin.jar'

def cross_with_catalogue(dir_path, cat_path):
    sex_files = [file.path for file in os.scandir(dir_path) if file.name.endswith('_sex.csv')]
    bat_path = os.path.join(dir_path, 'bat.ajs')
    for sex_path in sorted(sex_files):
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
        #os.system('java -jar C:\Aladin.jar -script < {}'.format(bat_path))


        filt = re.findall(r'.*_(\w+)_sex.csv$', sex_path.split(os.sep)[-1])[0].upper()

        try:
            os.remove(bat_path)

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
