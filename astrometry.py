import os
import re
from pyfits import open as pyfits_open


USE_SEXTRACTOR = '--use-sextractor '


def is_float_digit(value):
    """
    Checks whether the value is a floating-point number.
    :return: True if value is float else False
    """
    return re.findall(r'^[-+]?\d+\.?\d*$', value)


def astrometry(dir_path):
    """
    Runs the solve-field soft for the list of fit-files. If header of fit-file has not empty values of 'RA' and 'DEC',
    then to the solve-field will be uses these values. Renames new-files to fits-files (astrometry returns
    fits files with name: 'files.new'). Deletes excess files (astrometry returns many excess files).
    :param dir_path: path to directory with fit-files
    :return: fits-files with astrometrical calibrations located in dir_path folder
    """ 
    fit_files = [file for file in os.scandir(dir_path) if file.name.endswith('.fit')]
    for file in sorted(fit_files, key=lambda x: x.name):
        pre_filename = file.path.split('.')[0]
        try:
            with pyfits_open(file.path) as f:
                header = f[0].header
            if is_float_digit(header['RA']) and is_float_digit(header['DEC']):
                os.system('./solve-field --ra {} --dec {} --radius 0.4 --cpulimit 5 {}{}'.format(header['RA'],
                                                                                                 header['DEC'],
                                                                                                 USE_SEXTRACTOR,
                                                                                                 file.path))
            else:
                os.system('./solve-field {}{}'.format(USE_SEXTRACTOR, file.path))
        except PermissionError as e:
            print(e)
        except KeyError as e:
            print(e)
            os.system('./solve-field {}{}'.format(USE_SEXTRACTOR, file.path))
        except OSError as e:
            print(e)
        else:
            try:
                os.rename('{}.{}'.format(pre_filename, 'new'),
                          '{}.{}'.format(pre_filename, 'fits'))
                os.system('chmod 777 {}'.format('{}.{}'.format(pre_filename, 'fits')))
            except OSError as e:
                print(e)
            else:
                try:
                    os.remove(file.path)
                except OSError as e:
                    print(e)

        del_files = [del_file for del_file in os.scandir(dir_path)
                                           if del_file.is_file()
                                           and re.findall('{}'.format(pre_filename.split(os.sep)[-1]), del_file.name)
                                           and not del_file.name.endswith('.fit')
                                           and not del_file.name.endswith('.fits')]

        for del_file in del_files:
            try:
                os.remove(del_file.path)
            except OSError as e:
                print(e)


def start_astrometry():
    '''
    Asks the user for the address with the fit-files.Starts astrometry(dir_path, fit_files). 
    Changes the permissions of the folder with the fit-files.
    :return: results of astrometry(dir_path, fit_files) 
    '''
    dir_path = os.path.abspath(input('Enter the path to directory:\n'))
    astrometry(dir_path)
    os.system('chmod -R 777 {}'.format(dir_path))


if __name__ == '__main__':
    start_astrometry()
