Help on module nomus:

NAME
    nomus

FUNCTIONS
    add_all_stars(stars_XMatch, header_NOMAD1, stars_NOMAD1, header_USNO, stars_USNO)
        Adds stars to the XMatch catalogue for which no intersections have been found in the NOMAD1 and USNO-B1 catalogs.
        :param stars_XMatch: data about stars from XMatch catalogue
        :param header_NOMAD1: header of data from NOMAD1 catalogue
        :param stars_NOMAD1: data about stars from NOMAD1 catalogue
        :param header_USNO: header of data from USNO-B1 catalogue
        :param stars_USNO: data about stars from USNO-B1 catalogue
        :return:

    calc_today_coordinates(stars_from_file, epoch)
        Calculates the objects' coordinates, which are in the XMatch_NOMAD_USNO_path file for the epoch,
        indicated in the variable epoch.
        :param stars_from_file: XMatch NOMAD1 and USNO-B1 catalogues
        :param epoch: epoch
        :return: XMatch NOMAD1 and USNO-B1 catalogues with added colums of new coordinates

    column_from_list_dct(data, dct_key)
        Gets a list of data from the list of dictionaries by key in dictionaries.
        :param data: list of dictionaries
        :param dct_key: key in dictionaries
        :return: list of data on the key in dictionaries.

    convert_XMatch_to_catalogue(stars_inform, epoch='')
        Converts XMatch-file from Aladin to csv- catalogue.
        :param NOMAD_USNO_path: path to XMatch NOMAD1 and USNO-B1 catalogues
        :return: csv- catalogue

    get_XMatch_NOMAD_USNO(save_dir_path, ra, de, radius="14'")
        Gets XMatch-file from Aladin.
        :param save_dir_path: path to save directory for catalogue from Aladin
        :param ra: RA
        :param de: DE
        :param radius: radius of area for catalogue
        :return: path to XMatch-file from Aladin

    open_csv(csv_file_path, delimiter=',', first_line=False)
        Open csv file
        :param csv_file_path: Path to csv-file
        :param delimiter: delimiter in csv-file
        :param first_line: Skip (True) the first line or not (False)
        :return: header and data from csv-file

    replace_data_in_dct(dct, first_symbol, second_symbol)
        Replaces first_symbol with second_symbol in dictionary.
        :param dct: data from dictionary
        :param first_symbol: first symbol
        :param second_symbol: second symbol
        :return: data with replaced symbols

    start_nomus()
        Starts the program.
        :return:

    write_to_csv(save_path, data, header, delimiter=',')
        Writes data in csv-file.
        :param save_path: savefile path
        :param data: data
        :param header: header of data
        :param delimiter: delimiter
        :return: csv-file with data

DATA
    ALADIN_PATH = r'C:\Aladin.jar'

FILE
    https://github.com/Vasilii-Moskvin/Astronomy/blob/master/nomus/nomus.py

EXAMPLES
    >python nomus.py
    Enter the path where to save NomUs catalogue:
    C:\Users\vasil\Desktop\
    Enter Epoch (by default is 2000):
    2018
    Enter RA:
    00:04:00.21
    Enter DE:
    +00:59:56.3
    Enter radius (arc minute, by default is 14):
    30
    File C:\Users\vasil\Desktop\nomus_cat.csv has created!

SETTINGS
    Set path to Aladin.jar. 7 line in the file, the value of ALADIN_PATH.
    The input format for RA and DE is the same as in the Aladin's scripts.
