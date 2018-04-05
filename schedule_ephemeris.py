from collections import namedtuple
import csv
import os.path

import JDN


__author__ = 'Vasilii Moskvin https://github.com/Vasilii-Moskvin'

Ephemeris = namedtuple('Ephemeris', ('object_name', 'start', 'center', 'end'))


class VarStar:
    """
    Class VarStar (Variable Star) is designed to compile an ephemeris schedule.
    """
    def __init__(self, object_name, ra, de, period, center_epoch, duration):
        self.object_name = object_name
        self.ra = ra
        self.de = de
        self.period = float(period)
        self.center_epoch = float(center_epoch)
        self.duration = float(duration)

    def start_epoch(self, center_epoch):
        """
        Returns start epoch for the center epoch
        :param center_epoch: a center epoch
        :return: start epoch for the center epoch
        """
        return center_epoch - self.duration / 2

    def end_epoch(self, center_epoch):
        """
        Returns end epoch for the center epoch
        :param center_epoch: a center epoch
        :return: end epoch for the center epoch
        """
        return center_epoch + self.duration / 2

    def next_ephemeris(self, epoch):
        """
        Returns next ephemeris from the epoch
        :param epoch: epoch of moment
        :return: next ephemeris from the epoch
        """
        return epoch + self.period

    def ephemerises_between_dates(self, first_date, second_date):
        """
        Returns list of ephemerises between dates for the star
        :param first_date: first date
        :param second_date: second date
        :return: list of ephemerises between dates for the star
        """
        first_date = float(first_date)
        second_date = float(second_date)
        lst_ephemerises = list()

        count_periods = (first_date - self.start_epoch(self.center_epoch)) // self.period
        temp_center_epoch = self.center_epoch + self.period * count_periods
        while self.next_ephemeris(self.end_epoch(temp_center_epoch)) <= second_date:
            lst_ephemerises.append(Ephemeris(self.object_name,
                                             self.start_epoch(temp_center_epoch),
                                             temp_center_epoch,
                                             self.end_epoch(temp_center_epoch)))
            temp_center_epoch = self.next_ephemeris(temp_center_epoch)

        return lst_ephemerises

    def __str__(self):
        return 'name = {} ra = {} de = {} period = {} center_epoch = {} duration = {}'.format(self.name, self.ra,
                                                                                              self.de, self.period,
                                                                                              self.center_epoch,
                                                                                              self.duration)


def stars_from_file(file_path):
    """
    Returns stars from file.
    File format - csv. Field names: object_name,ra,de,period,center_epoch,duration

    Example of file:

    object_name,ra,de,period,center_epoch,duration
    HAT-P12b,13 57 33.684,+43 29 37.35,3.2130598,2457773.6298611113,0.09743
    WASP-12b,06 30 32.79,+29 40 20.4,1.0914222,2457773.419444444,0.12504

    :param file_path: path to file
    :return: stars from file
    """
    with open(file_path, 'r') as csv_file:
        csv_file = csv.DictReader(csv_file, delimiter=',')
        lst_stars = [VarStar(**row) for row in csv_file]

    return lst_stars


def ephemerises_for_stars(stars, first_date, second_date):
    """
    Returns sorted list of ephemerises for star between first_date and second_date
    :param stars: list of stars
    :param first_date: first date
    :param second_date: second date
    :return: sorted list of ephemerises for star between first_date and second_date
    """
    lst_ephemerises = list()
    for star in stars:
        lst_ephemerises.extend(star.ephemerises_between_dates(first_date, second_date))

    return sorted(lst_ephemerises, key=lambda x: x.start)


def print_ephemerises(lst_ephemerises):
    """
    Print the list of ephemerises
    :param lst_ephemerises: list of ephemerises
    :return: printed list of ephemerises
    """
    previous_JD = lst_ephemerises[0].start
    previous_hour = JDN.get_GD(str(lst_ephemerises[0].start))[3]
    previous_day = JDN.get_GD(str(lst_ephemerises[0].start))[2]
    print('{}/{:02}/{:02}'.format(*JDN.get_GD(str(lst_ephemerises[0].start))))
    for ephemeris in lst_ephemerises:
        this_JD = ephemeris.start
        this_hour = JDN.get_GD(str(ephemeris.start))[3]
        this_day = JDN.get_GD(str(ephemeris.start))[2]
        if this_JD - previous_JD >= 0.5:
            if this_hour >= 12:
                print('\n{}/{:02}/{:02}'.format(*JDN.get_GD(str(ephemeris.start))))
                previous_JD = this_JD
                previous_hour = this_hour
                previous_day = this_day
            else:
                if this_day != previous_day:
                    print('\n{}/{:02}/{:02}'.format(*JDN.get_GD(str(ephemeris.start - 1))))
                    previous_JD = this_JD
                    previous_hour = this_hour
                    previous_day = this_day
        else:
            if this_hour >= 12:
                if previous_hour < 12:
                    print('\n{}/{:02}/{:02}'.format(*JDN.get_GD(str(ephemeris.start))))
                    previous_JD = this_JD
                    previous_hour = this_hour
                    previous_day = this_day
        print("{:10}\t"
              "{}/{:02}/{:02} {:02}:{:02}:{:02}\t"
              "{}/{:02}/{:02} {:02}:{:02}:{:02}\t"
              "{}/{:02}/{:02} {:02}:{:02}:{:02}".format(ephemeris.object_name,
                                                        *JDN.get_GD(str(ephemeris.start)),
                                                        *JDN.get_GD(str(ephemeris.center)),
                                                        *JDN.get_GD(str(ephemeris.end))))


def main():
    stars = stars_from_file(os.path.abspath(input('Enter path:\n')))
    first_date = JDN.get_JD(*input('Enter first date'
                                   ' (format: year/month/day/hour/minute/second):\n').strip().split('/'))
    second_date = JDN.get_JD(*input('Enter first date '
                                    '(format: year/month/day/hour/minute/second):\n').strip().split('/'))
    lst_ephemerises = ephemerises_for_stars(stars, str(first_date), str(second_date))
    print_ephemerises(lst_ephemerises)


if __name__ == '__main__':
    main()
