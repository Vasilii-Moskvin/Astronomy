import os
import re
from func.work_with_csv import write_csv


def sextractor(fits_files):
    """
    Runs s-extractor for files from the dirpath
    :param fits_files: list of information about fits-files
    :return: the result of the s-extractor
    """
    for file in sorted(fits_files, key=lambda x: x.name):
        file_prename = re.findall(r'^(.+)\.fits$', file.path)[0]
        sex_path = '{}_sex.txt'.format(file_prename)
        try:
            os.system('sextractor -CATALOG_NAME {} -c daofind.sex {}'.format(sex_path, file.path))
        except OSError as e:
            print(e)
        except UnicodeDecodeError as e:
            print(e)

        data = []
        units = []
        pre_filename = re.findall(r'^(.+_sex)\.txt$', sex_path)[0]

        try:
            with open(sex_path, 'r') as f:
                for line in f:
                    if line[0] == '#':
                        units.append(line.strip().split()[2])
                    else:
                        data.append(line.strip().split())
            
            write_csv(sex_path, units, data)
            os.rename('{}.{}'.format(pre_filename, 'txt'),
                      '{}.{}'.format(pre_filename, 'csv'))
            os.system('chmod 777 {}'.format('{}.{}'.format(pre_filename, 'csv')))
        except OSError as e:
            print(e)


def start_sextractor():
    '''
    Asks the user for the address with the fits-files. Starts sextractor(dir_path). 
    Changes the permissions of the folder with the fits-files.
    :return: catalogues of objects from fits-files plased in the dir_path.
    '''
    dir_path = os.path.abspath(input('Enter the path to directory:\n'))
    fits_files = [file for file in os.scandir(dir_path) if file.name.endswith('.fits')]
    sextractor(fits_files)
    os.system('chmod -R 777 {}'.format(dir_path))


if __name__ == '__main__':
    start_sextractor()
