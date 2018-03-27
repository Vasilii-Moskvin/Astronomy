Help on module get_catalogue:

NAME
    get_catalogue

FUNCTIONS
    get_catalogue()
        Reades fit-files and *._cat.csv in the folder dir_path. Produces astrometric calibration over the fit-files.
        Produces photometry over the fits-files in the SExtractor software. Produces a compilation catalogue of objects
        from data received after processing in the SExtractor software. Produces a correlation analysis of data in this catalogue,
        to find and take into account distortions obtained in the process of astronomical observations. Reduction of the received data
        to the reference catalog system, which is in the file *._cat.csv.
        :return: istrumental and reduced catalogs, reduction lines, correlograms for a particular filter.

FILE
    https://github.com/Vasilii-Moskvin/Astronomy/blob/master/get_catalogue/get_catalogue.py

EXAMPLES
    >python get_catalogue.py
    Enter the path to directory:
    /home/crao/Astronomy/path_to_fit_files
    
    >python get_catalogue.py
    Enter the path to directory:
    /home/crao/Astronomy/path_to_fit_files --two

SETTINGS:
    Copy get_catalogue.py and func/ to folder with solve-field file (nova.astrometry.net software).
    Copy the fit-files and the reference catalogue to the folder path_to_fit_files. 
    Run the script get_catalogue.py. Specify the path to the directory with the fit-files
    and the reference directory. If the flag --two is not specified, then only those objects 
    that are in the reference catalogue are searched. If this flag is specified,
    then search for all objects that are in at least two frames.

    Format of the name of the fit-file:
    WaSP12-001_R.fit
    WASP12_01-011_V.fit

    Format of the reference catalohue:
    Format of the name of the header:
    RAJ2000,DEJ2000,NOMAD1_ID,USNO-B1.0,B_mag,V_mag,R_mag,I_mag
    Format of the empty values:
    -
    Format of the coodinates:
    Like format of the coodinates in the Aladin software.
    If you use the --two flag, it is preferable to use coordinates in the floating-point number format.
    