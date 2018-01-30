from collections import OrderedDict, namedtuple
import csv
import os.path


def get_XMatch_NOMAD_USNO(save_dir_path, ra, de, radius="14'"):
    """
    Gets XMatch-file from Aladin.
    :param save_dir_path: path to save directory for catalogue from Aladin
    :param ra: RA
    :param de: DE
    :param radius: radius of area for catalogue
    :return: path to XMatch-file from Aladin
    """
    XMatch_NOMAD_USNO_path = os.path.join(save_dir_path, 'nomus_cat.csv')
    script_for_Aladin = ['#AJS',
                         'sync',
                         'get VizieR(NOMAD1) {} {} {}'.format(ra, de, radius),
                         'get VizieR(USNO-B1) {} {} {}'.format(ra, de, radius),
                         'xmatch NOMAD1 USNO-B1 1.0',
                         'export XMatch {}'.format(XMatch_NOMAD_USNO_path),
                         'quit']

    script_for_Aladin_path = os.path.join(save_dir_path, 'script_for_Aladin.txt')
    with open(script_for_Aladin_path, 'w') as w:
        for line in script_for_Aladin:
            w.write('{}\n'.format(line))

    os.system('java -jar C:\Aladin.jar -script < {}'.format(script_for_Aladin_path))
    os.remove(script_for_Aladin_path)

    if os.path.exists(XMatch_NOMAD_USNO_path):
        return XMatch_NOMAD_USNO_path
    else:
        return False


def convert_XMatch_to_catalogue(XMatch_NOMAD_USNO_path):
    """
    Converts XMatch-file from Aladin to csv- catalogue.
    :param NOMAD_USNO_path: path to XMatch NOMAD1 and USNO-B1 catalogues
    :return: csv- catalogue
    """
    convert_header = OrderedDict((('_RAJ2000_tab1', 'RAJ2000'),
                                  ('_DEJ2000_tab1', 'DEJ2000'),
                                  ('Bmag_tab1', 'B_mag'),
                                  ('Vmag_tab1', 'V_mag'),
                                  ('Rmag_tab1', 'R_mag'),
                                  ('Imag_tab2', 'I_mag')))
    Stars = namedtuple('Stars', convert_header.values())

    with open(XMatch_NOMAD_USNO_path, 'r') as csv_file:
        csv_file = csv.DictReader(csv_file, delimiter='\t')
        next(csv_file)
        stars_from_file = [Stars(**dict(map(lambda key: (convert_header[key], star[key]) if star[key] != ' ' else
                                                        (convert_header[key], '-'),
                                            filter(lambda x: x in convert_header, star.keys()))))
                           for star in csv_file]

    with open(XMatch_NOMAD_USNO_path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        writer.writerow(convert_header.values())
        writer.writerows(stars_from_file)


def calc_today_coordinates(NOMAD_USNO_path, epoch):
    pass


def main():
    save_dir_path = os.path.abspath(input('Enter the path where to save NomUs catalogue:\n'))
    epoch = input('Enter Epoph:\n')
    ra = input('Enter RA:\n')
    de = input('Enter DE:\n')
    radius = input('Enter radius:\n')

    try:
        XMatch_NOMAD_USNO_path = get_XMatch_NOMAD_USNO(save_dir_path, ra, de, radius)
        if XMatch_NOMAD_USNO_path:
            convert_XMatch_to_catalogue(XMatch_NOMAD_USNO_path)
    except OSError as e:
        print(e)


if __name__ == '__main__':
    main()
