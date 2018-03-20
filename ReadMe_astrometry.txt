Help on module astrometry:

NAME
    astrometry

FUNCTIONS
    astrometry(dir_path, fit_files)
        Runs the solve-field soft for the list of fit-files. If header of fit-file has not empty values of 'RA' and 'DEC',
        then to the solve-field will be uses these values. Renames new-files to fits-files (astrometry returns
        fits files with name: 'files.new'). Deletes excess files (astrometry returns many excess files).
        :param dir_path: path to directory with fit_files
        :param fit_files: path to fit-files
        :return: fits-files with astrometrical calibrations located in dir_path folder

    is_float_digit(value)
        Checks whether the value is a floating-point number.
        :return: True if value is float else False

    start_astrometry()
        Asks the user for the address with the fit-files.Starts astrometry(dir_path, fit_files).
        Changes the permissions of the folder with the fit-files.
        :return: results of astrometry(dir_path, fit_files)

DATA
    USE_SEXTRACTOR = '--use-sextractor '

FILE
    https://github.com/Vasilii-Moskvin/Astronomy/blob/master/astrometry.py

EXAMPLES
    >python astrometry.py
    Enter the path to folder with fit-files:
    /home/crao/Astronomy/Processing/test

SETTINGS:
    Copy astrometry.py to folder with solve-field file.
    If you want to use SExtractor in the work of solve-field set: USE_SEXTRACTOR = '--use-sextractor '.
    If you don't want to use SExtractor in the work of solve-field set: USE_SEXTRACTOR = ''.
