import os.path
from func.calc_astrometry import astrometry
from func.JD_form_fits_files import JD_form_fits_files
from func.run_sextractor import sextractor
from func.cross_with_catalogue import cross_with_catalogue
from func.full_inst_cat import full_inst_cat
from func.build_catalogue import load_catalogue, build_catalogues, write_catalogue
from func.corr_analysis import produce_corr_analysis
from func.make_reduction import make_reducted_catalogue, make_reducted_catalogue_1, get_ABC, save_regress_line
from func.graphics import save_correlogram, save_regress
from func.StarFromCSV import StarFromCSV
from func.work_with_csv import open_dict_csv, write_dict_csv


def get_catalogue():
    '''
    Reades fit-files and *._cat.csv in the folder dir_path. Produces astrometric calibration over the fit-files.
    Produces photometry over the fits-files in the SExtractor software. Produces a compilation catalogue of objects
    from data received after processing in the SExtractor software. Produces a correlation analysis of data in this catalogue,
    to find and take into account distortions obtained in the process of astronomical observations. Reduction of the received data
    to the reference catalog system, which is in the file *._cat.csv.
    :return: istrumental and reduced catalogs, reduction lines, correlograms for a particular filter.
    '''
    input_txt = input('Enter the path to directory:\n').strip()
    input_lst = input_txt.split()
    two_flag = False
    if len(input_lst) == 2 and input_lst[1] == '--two':
        two_flag = True
    dir_path = os.path.abspath(input_lst[0])
    #obj = input('Enter the RA and DE of object (RA_DE):\n')
    obj = '303.381794_+65.162131'
    
    fit_files = [file for file in os.scandir(dir_path) if file.name.endswith('.fit')]
    fits_files = astrometry(dir_path, fit_files)
    JD_path = JD_form_fits_files(dir_path, fits_files)
    sextractor(fits_files)
    cat_path = [file.path for file in os.scandir(dir_path) if file.name.endswith('_cat.csv')][0]
    cross_with_catalogue(dir_path, cat_path, two_flag)
    full_inst_cat_path = full_inst_cat(dir_path, JD_path)

    data, prefix_full_cat = load_catalogue(dir_path, full_inst_cat_path)
    to_sys_data = build_catalogues(data, 'to_sys')
    # write_catalogue(dir_path, to_sys_data, prefix_full_cat, 'inst', 'to_sys')

    data_corr, to_sys_data_corr = produce_corr_analysis(dir_path, data, to_sys_data, obj)
    write_catalogue(dir_path, data_corr, prefix_full_cat, 'corr', 'full')
    # write_catalogue(dir_path, to_sys_data_corr, prefix_full_cat, 'corr', 'to_sys')

    if len(StarFromCSV.filts) == 1:
        reducted_catalogue, ABC_by_filt = make_reducted_catalogue(to_sys_data_corr)
        save_regress_line(dir_path, reducted_catalogue, ABC_by_filt)
        full_reduction_cat = make_reducted_catalogue_1(data_corr, ABC_by_filt)
        write_catalogue(dir_path, full_reduction_cat, prefix_full_cat, '', 'result')
        two_reduction_cat = make_reducted_catalogue_1(data_corr, ABC_by_filt, two=True)
        write_catalogue(dir_path, two_reduction_cat, prefix_full_cat, 'two', 'result')
    #else:
    #    write_catalogue(dir_path, reducted_catalogue, prefix_full_cat, '', 'result')


    StarFromCSV.reset_class_data()


if __name__ == '__main__':
    get_catalogue()
