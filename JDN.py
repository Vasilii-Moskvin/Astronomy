from math import modf

__author__ = "Moskvin Vasilii (vasiliy.moscvin@yandex.ru)"

def get_leap(year):
    """
    Returns True if year is Leap year, else False
    :param year: year
    :return: True if year is Leap year, else False
    """
    if year % 400 == 0:
        return True
    elif year % 100 == 0:
        return False
    elif year % 4 == 0:
        return True
    else:
        return False

def check_output_data(year, month):
    """
    Returns the number of days in the month
    :param year: year
    :param month: month
    :return: days in the month
    """
    dct_days_in_months = {True: {1: 31,
                                 2: 29,
                                 3: 31,
                                 4: 30,
                                 5: 31,
                                 6: 30,
                                 7: 31,
                                 8: 31,
                                 9: 30,
                                 10: 31,
                                 11: 30,
                                 12: 31},
                           False: {1: 31,
                                   2: 28,
                                   3: 31,
                                   4: 30,
                                   5: 31,
                                   6: 30,
                                   7: 31,
                                   8: 31,
                                   9: 30,
                                   10: 31,
                                   11: 30,
                                   12: 31}}
    return dct_days_in_months[get_leap(year)][month]

def get_JDN(year, month, day):
    """
    Return  Julian Day Number (JDN)
    :param year: year
    :param month: month
    :param day: day
    :return: JDN
    """
    a = (14 - int(month)) // 12
    y = int(year) + 4800 - a
    m = int(month) + 12 * a - 3
    JDN = int(day) + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    return JDN


def get_JD(year, month, day, h, m, s):
    """
    Return Julian date (JD)
    :param year: year
    :param month: month
    :param day: day
    :param h: hour
    :param m: minute
    :param s: second
    :return: JD
    """
    return get_JDN(year, month, day) + (float(h) - 12) / float(24) + float(m) / 1440 + float(s) / 86400


def get_GD(JD, *key):
    """
    Return Gregory date (GD) from Julian date (JD)
    :param JD: Julian date
    :param key: If key == -f it return formated Gregory date
    :return: Gregory date
    """
    Z = int(JD.split(".")[0])
    F = round(float(JD) - Z, 6)
    if Z < 2299161:
        A = Z
    else:
        alfa = (Z - 1867216.25) // 36524.25
        A = Z + 1 + alfa - alfa // 4
    B = A + 1524
    C = (B - 122.1) // 365.25
    D = int(365.25 * C)
    E = (B - D) // 30.6001
    ost, day = modf(B - D - int(30.6001 * E) + F)
    day = int(day)
    hours = int(ost * 24)
    minuts = int((ost - hours / 24) * 1440)
    seconds = int((ost - hours / 24 - minuts / 1440) * 86400)
    if E < 13.5:
        month = int(E - 1)
    else:
        month = int(E - 13)
    if month > 2.5:
        year = int(C - 4716)
    else:
        year = int(C - 4715)
    if hours > 11:
        if day == check_output_data(year, month):
            if month == 12:
                month = 1
                year += 1
            else:
                month += 1
            day = 1
        else:
            day += 1
        hours -= 12
    else:
        hours += 12
    if key and key[0] == '-f':
        return "{}/{}/{}\t{}:{}:{}".format(day, month, year, hours, minuts, seconds)
    else:
        return year, month, day, hours, minuts, seconds


if __name__ == '__main__':
    JD = input('Enter JD:\n')
    print(get_GD(JD, '-f'))
