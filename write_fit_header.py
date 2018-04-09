from astropy.io import fits
import os.path
import json


def load_inf_to_header(path_to_inf_to_header):
    """
    Loads information for fits header from path_to_inf_to_header
    :param path_to_inf_to_header: path to file with information for fits header
    :return: information for fits header
    """
    with open(path_to_inf_to_header, encoding='utf-8') as json_file:
        inf_to_header = json.load(json_file)

    return inf_to_header


def write_inf_to_header(header, inf_to_header):
    """
    Writes information to fits header
    :param header: header of fits file
    :param inf_to_header: information to fits header
    :return: changed header variable
    """
    for key, value in inf_to_header.items():
        header[key] = value


def change_fits_files(dir_fits_path, inf_to_header):
    """
    Changes fits header. Gets list of fit\fits files in dir_fits_path and re-writes fits header
    :param dir_fits_path: path to dir with fit\fits files
    :param inf_to_header: information to fits header
    :return: changed fit\fits files
    """
    lst_fits_files = [file.path for file in os.scandir(dir_fits_path)
                      if file.name.endswith('.fits') or file.name.endswith('.fit')]

    for fits_file_path in lst_fits_files:
        data, header = fits.getdata(fits_file_path, header=True)
        for key, value in inf_to_header.items():
            header[key] = value

        fits.writeto(fits_file_path, data, header=header, overwrite=True)
        print('{} Ok!'.format(fits_file_path.split()[-1]))


def main():
    dir_fits_path = os.path.abspath(input('Enter path to directory with fits-files:\n'))
    path_to_inf_to_header = os.path.abspath(input('Enter path to inf to header data:\n'))

    inf_to_header = load_inf_to_header(path_to_inf_to_header)
    change_fits_files(dir_fits_path, inf_to_header)


if __name__ == '__main__':
    main()