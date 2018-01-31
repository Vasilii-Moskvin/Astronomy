from collections import OrderedDict, namedtuple
import csv
import os.path
import math
import often_used as ou


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

    os.system('java -jar C:\Aladin.jar -script < {}'.format(script_for_Aladin_path))
    os.remove(script_for_Aladin_path)

    if os.path.exists(XMatch_NOMAD_USNO_path):
        return XMatch_NOMAD_USNO_path
    else:
        return False


def convert_XMatch_to_catalogue( stars_inform, epoch=''):
    """
    Converts XMatch-file from Aladin to csv- catalogue.
    :param NOMAD_USNO_path: path to XMatch NOMAD1 and USNO-B1 catalogues
    :return: csv- catalogue
    """
    if epoch:
        convert_header = OrderedDict((('_RAJ{}_tab1'.format(epoch), 'RAJ{}'.format(epoch)),
                                      ('_DEJ{}_tab1'.format(epoch), 'DEJ{}'.format(epoch)),
                                      ('NOMAD1_tab1', 'NOMAD1_ID'),
                                      ('USNO-B1.0_tab2', 'USNO-B1.0'),
                                      ('Bmag_tab1', 'B_mag'),
                                      ('Vmag_tab1', 'V_mag'),
                                      ('Rmag_tab1', 'R_mag'),
                                      ('Imag_tab2', 'I_mag')))
    else:
        convert_header = OrderedDict((('_RAJ2000_tab1', 'RAJ2000'),
                                      ('_DEJ2000_tab1', 'DEJ2000'),
                                      ('NOMAD1_tab1', 'NOMAD1_ID'),
                                      ('USNO-B1.0_tab2', 'USNO-B1.0'),
                                      ('Bmag_tab1', 'B_mag'),
                                      ('Vmag_tab1', 'V_mag'),
                                      ('Rmag_tab1', 'R_mag'),
                                      ('Imag_tab2', 'I_mag')))

    save_header = list(convert_header.values())

    data_to_csv = [dict(map(lambda key: (convert_header[key], star[key]),
                            filter(lambda x: x in convert_header, star.keys()))) for star in stars_inform]

    return data_to_csv, save_header


def calc_today_coordinates(stars_from_file, epoch):
    '''
    Calculates the objects' coordinates, which are in the XMatch_NOMAD_USNO_path file for the epoch,
    indicated in the variable epoch.
    :param stars_from_file: XMatch NOMAD1 and USNO-B1 catalogues
    :param epoch: epoch
    :return: XMatch NOMAD1 and USNO-B1 catalogues with added colums of new coordinates
    '''
    ra = '_RAJ2000_tab1'
    dec = '_DEJ2000_tab1'
    pm_ra = 'pmRA_tab1'
    pm_dec = 'pmDE_tab1'
    ra_by_epoch = '_RAJ{}_tab1'.format(epoch)
    dec_by_epoch = '_DEJ{}_tab1'.format(epoch)

    float_epoch = float(epoch) - 2000

    for star in stars_from_file:
        star[dec_by_epoch] = float(star[dec]) + (0.001 * float(star[pm_dec]) * float_epoch) / 3600
        old_cos_dec = math.cos((math.pi / 180) * float(star[dec]))
        new_cos_dec = math.cos((math.pi / 180) * float(star[dec_by_epoch]))
        star[ra_by_epoch] = float(star[ra]) + \
                            ((0.001 * float(star[pm_ra]) * float_epoch) / 3600) * (new_cos_dec / old_cos_dec)

    return stars_from_file


def main():
    save_dir_path = os.path.abspath(input('Enter the path where to save NomUs catalogue:\n'))
    epoch = input('Enter Epoch:\n')
    ra = input('Enter RA:\n')
    de = input('Enter DE:\n')
    radius = input('Enter radius:\n')

    try:
        XMatch_NOMAD_USNO_path = get_XMatch_NOMAD_USNO(save_dir_path, ra, de, radius)
        if XMatch_NOMAD_USNO_path:
            header, stars_inform = ou.open_csv(XMatch_NOMAD_USNO_path, delimiter='\t', first_line=True)
            stars_inform = ou.replace_data_in_dcts(stars_inform, ' ', '-')
            if epoch:
                stars_inform = calc_today_coordinates(stars_inform, epoch)
            stars_inform, save_header = convert_XMatch_to_catalogue(stars_inform, epoch)
            ou.write_to_csv(XMatch_NOMAD_USNO_path, stars_inform, save_header, delimiter=',')
    except OSError as e:
        print(e)


if __name__ == '__main__':
    main()
