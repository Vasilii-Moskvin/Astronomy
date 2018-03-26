import copy
import re
import numpy as np
from collections import OrderedDict


class StarFromCSV(OrderedDict):
    '''
    Works with information about the stars in the compilation directory. Work with data can also be done as OrderedDict.
    Data in StarFromCSV is divided into filters, in which frames were made. In each filter, the data is divided into two parts.
    In the first series of stellar magnitudes. The second statistics of these series. In this statistical part are the minimum,
    maximum, average stellar magnitude, error of stellar magnitudes, the number of frames on which the object is observed. Data
    in the statistical part calculating by first part.
    '''
    header = []
    filts = []
    n_fields = OrderedDict()
    mag_fields_name = OrderedDict()
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.get_static_info(*args)

    @classmethod
    def reset_class_data(cls):
        '''
        Resets values for class variables StarFromCSV.
        :return: reseted values: header, filts, n_fields, mag_fields_name.
        '''
        cls.header = []
        cls.filts = []
        cls.n_fields = OrderedDict()
        cls.mag_fields_name = OrderedDict()

    @classmethod
    def get_static_info(cls, *args):
        '''
        Fills class information about filters, frame numbers, and statistical fields.
        :return: filled class information about filters, frame numbers, and statistical fields.
        '''
        if not cls.header:
            cls.header = [item[0] for item in args[0]]
        if not cls.filts:
            cls.filts = list(map(lambda x: re.findall(r'(^\w+)_n$', x)[0],
                              filter(lambda x: re.findall(r'^\w+_n$', x), cls.header)))
        if not cls.n_fields:
            cls.n_fields = OrderedDict(map(lambda field: (field, ''),
                                           filter(lambda x: re.findall(r'^\w+_n$', x), cls.header)))
        if not cls.mag_fields_name:
            for filt in cls.filts:
                cls.mag_fields_name[filt] = list(filter(lambda x: re.findall(r'{}_mag_\d+$'.format(filt), x),
                                                        cls.header))

    @classmethod
    def fill_n_fields(cls, dct_max_n_value):
        '''
        Updates maximum number of frames on which a object is located is determined in the all filters.
        :param dct_max_n_value: the dictionary of maximum number of frames on which a object is located
                                is determined in the all filters.
        :return: updated maximum number of frames on which a object is located is determined in the all filters.
        '''
        for field in cls.n_fields.keys():
            cls.n_fields[field] = dct_max_n_value[field]

    @property
    def mag_values(self):
        '''
        Returns data with serieses of stellar magnitudes, depending on the filter.
        :return: data with serieses of stellar magnitudes for a particular filter.
        '''
        mag_values = OrderedDict([(filt, OrderedDict([(key, value) for key, value in self.items()
                                                      if re.findall(r'{}_mag_\d+$'.format(filt), key)]))
                                  for filt in StarFromCSV.filts])
        return mag_values

    def mag_digit_values(self, filt):
        '''
        Returns non-empty data with series of stellar magnitudes for a particular filter.
        :param filt: a particular filter.
        :return: non-empty data with series of stellar magnitudes for a particular filter.
        '''
        return filter(lambda x: x != '-', tuple(self.mag_values[filt].values()))

    def reduction_AC(self, ABC_by_filt):
        '''
        Produces a reduction of instrumental stellar magnitudes to the system of the selected catalog.
        Data on the reduction coefficients are in the dictionary ABC_by_filt.
        :param ABC_by_filt: the dictionary with the reduction coefficients.
        :return: the catalogue with a reduction of instrumental stellar magnitudes.
        '''
        filt = StarFromCSV.filts[0]
        a = ABC_by_filt[filt].a
        c = ABC_by_filt[filt].c

        for key, value in self.mag_values[filt].items():
            if value != '-':
                self[key] = a * value + c

        self.update_stat_mag_by_filter(filt)

    def reduction_ABC(self, ABC_by_filt):
        '''
        Produces a multi-color reduction of instrumental stellar magnitudes to the system of the selected catalog.
        Data on the reduction coefficients are in the dictionary ABC_by_filt.
        :param ABC_by_filt: the dictionary with the reduction coefficients for the all filters.
        :return: the catalogue with a multi-color reduction of instrumental stellar magnitudes.
        '''
        mag_values = copy.deepcopy(self.mag_values)
        for filt1 in StarFromCSV.filts:
            a = ABC_by_filt[filt1].a
            b = ABC_by_filt[filt1].b
            c = ABC_by_filt[filt1].c
            filt2 = self.choose_filter(filt1)
            for key_filt1, key_filt2 in zip(mag_values[filt1], mag_values[filt2]):
                if self[key_filt1] != '-' and self[key_filt1] != '-':
                    self[key_filt1] = a * mag_values[filt1][key_filt1] + b * (mag_values[filt2][key_filt2] - mag_values[filt1][key_filt1]) + c
                else:
                    self[key_filt1] = '-'

            self.update_stat_mag_by_filter(filt1)

    def mag_value_by_filter(self, filt):
        '''
        Returns data with series of stellar magnitudes for a particular filter.
        :param filt: a particular filter.
        :return: data with series of stellar magnitudes for a particular filter.
        '''
        try:
            output = self.mag_values[filt]
        except KeyError as e:
            print("Filter's name error. You input {} filter's name ".format(e))
        else:
            return output

    def change_mag(self, change_lst, filt):
        '''
        Changes the value of the star magnitudes to the values listed in the change_lst for a particular filter.
        :param change_lst: list with changes in stellar magnitudes.
        :param filt: a particular filter.
        :return: changed stellar magnitudes for a particular filter.
        '''
        for field, change in zip(self.mag_value_by_filter(filt).keys(), change_lst):
            if self[field] != '-':
                self[field] = float(self[field]) - change

        self.update_stat_mag_by_filter(filt)

    def update_stat_mag_by_filter(self, filt):
        '''
        Updates the statistical part of the object for a particular filter.
        :param filt: a particular filter.
        :return: Updated the statistical part of the object for a particular filter.
        '''
        mags = tuple(filter(lambda x: x != '-', self.mag_digit_values(filt)))
        if mags:
            self['{}_min_mag'.format(filt)] = max(mags)
            self['{}_max_mag'.format(filt)] = min(mags)
            self['{}_avr_mag'.format(filt)] = np.mean(mags)
            self['{}_std_mag'.format(filt)] = np.std(mags)
        else:
            self['{}_min_mag'.format(filt)] = '-'
            self['{}_max_mag'.format(filt)] = '-'
            self['{}_avr_mag'.format(filt)] = '-'
            self['{}_std_mag'.format(filt)] = '-'

    def update_mag(self):
        '''
        Updates the statistical part of the object in the all filters.
        :return: updated the statistical part of the object in the all filters.
        '''
        for filt in StarFromCSV.filts:
            self.update_stat_mag_by_filter(filt)

    @staticmethod
    def choose_filter(filt):
        '''
        Selects a filter to carry out multi-color reduction for a particular filter.
        :return: a filter to carry out multi-color reduction for a particular filter.
        '''
        if filt == 'V':
            filt2 = 'B'
        elif filt == 'R':
            filt2 = 'I'
        elif filt == 'B':
            filt2 = 'V'
        elif filt == 'I':
            filt2 = 'R'

        return filt2
