from func.graph_flare import draw_full_frame, draw_aperture, draw_LC, draw_LC_1, draw_LC_2
from func.Eat import Eat
from collections import namedtuple
import math
import csv
import numpy as np
import json
import os
import matplotlib.pyplot as plt
from pprint import pprint
from numba import jit
from PIL import Image
import JDN
import re
import logging
import imageio

logging.basicConfig(filename="logger.log", filemode='w', level=logging.DEBUG)


Photometry = namedtuple('Photometry', ('x', 'y', 'radius_star', 'gap', 'width_bkgnd'))
RAW_Curve = namedtuple('RAW_Curve', ('time', 'x', 'y'))
Curve = namedtuple('Curve', ('time', 'x', 'y', 'in_star'))
JD_data = namedtuple('JD_data', ('name', 'year', 'month', 'day', 'hour', 'min', 'sec'))

def set_photometry(xy, path_json):
    ans = ''
    while True:
        ans = input(r'Photometry settings (file\auto): ').strip()
        if ans == 'file':
            photometry_settings = Photometry(**json.load(open(path_json)))
            break
        elif ans == 'auto aperture':
            default_settings = Photometry(**json.load(open(path_json)))
            photometry_settings = calc_auto_aperture(xy, default_settings.x, default_settings.y)
            break
        elif ans == 'auto':
            photometry_settings = calc_total_auto_aperture(xy)
            break
        elif ans == 'q':
            photometry_settings = None
            break
        else:
            mprint("Wrong key specified. Enter file or auto")

    return photometry_settings


def get_std(xy):
    X = list(map(lambda x: sum(x), zip(*xy)))
    Y = list(map(lambda x: sum(x), xy))
    x_n = []
    for index, src in enumerate(X):
        for i in range(int(src)):
            x_n.append(index)
    y_n = []
    for index, src in enumerate(Y):
        for i in range(int(src)):
            y_n.append(index)
    std_x = np.std(x_n)
    std_y = np.std(y_n)

    return std_x, std_y


def calc_auto_aperture(xy, x0, y0):
    std_X, std_Y = get_std(xy)
    r_star = std_X if std_X > std_Y else std_Y
    gap = 0.44 * r_star
    width_bkgnd = 0.44 * r_star
    photometry_settings = Photometry(x0, y0, r_star, gap, width_bkgnd)

    return photometry_settings


def calc_total_auto_aperture(xy):
    y_0_star, x_0_star = tuple(map(int, np.unravel_index(xy.argmax(), xy.shape)))
    photometry_settings = calc_auto_aperture(xy, x_0_star, y_0_star)
    
    return photometry_settings


def circle_coord(photo):
    fi = [i for i in range(360)]
    x_star = [photo.radius_star * math.cos(i) + photo.x for i in fi]
    y_star = [photo.radius_star * math.sin(i) + photo.y for i in fi]
    x_background_inside = [(photo.radius_star + photo.gap) * math.cos(i) + photo.x for i in fi]
    y_background_insede = [(photo.radius_star + photo.gap) * math.sin(i) + photo.y for i in fi]
    x_background_outside = [(photo.radius_star + photo.gap + photo.width_bkgnd) * math.cos(i) + photo.x for i in fi]
    y_background_outsede = [(photo.radius_star + photo.gap + photo.width_bkgnd) * math.sin(i) + photo.y for i in fi]

    return x_star, y_star, x_background_inside, y_background_insede, x_background_outside, y_background_outsede


def write_dict_csv(save_path, header, data, delimiter=','):
    '''
    Writes data in csv-file.
    :param save_path: savefile path
    :param data: data
    :param header: header of data
    :param delimiter: delimiter
    :return: csv-file with data
    '''

    with open(save_path, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=header, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(data)
    if os.path.exists(save_path):
        mprint('File {} has created!'.format(save_path))
    else:
        mprint('Error: File {} has no created!'.format(save_path))


def open_dict_csv(csv_file_path, delimiter=',', first_line=False):
    '''
    Open csv file
    :param csv_file_path: Path to csv-file
    :param delimiter: delimiter in csv-file
    :param first_line: Skip (True) the first line or not (False)
    :return: header and data from csv-file
    '''
    with open(csv_file_path, 'r') as csv_file:
        csv_file = csv.DictReader(csv_file, delimiter=delimiter)
        header = csv_file.fieldnames
        if first_line:
            next(csv_file)
        data = [row for row in csv_file]

    return header, data


def save_json_photo(json_path, settings_photometry):
    with open(json_path, "w", encoding="utf-8") as file:
        json.dump(settings_photometry._asdict(), file)
    if os.path.exists(json_path):
        mprint('File {} has created!'.format(json_path))
    else:
        mprint('Error: File {} has no created!'.format(json_path))


def save_graph_lc(save_path, time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta):
    draw_LC(time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta)
    plt.savefig(save_path,
                  dpi=400,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None,
                  fmt='png')
    plt.close()
    if os.path.exists(save_path):
        mprint('File {} has created!'.format(save_path))
    else:
        mprint('Error: File {} has no created!'.format(save_path))

def save_graph_lc_1(save_path, JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta):
    draw_LC_1(JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta)
    plt.savefig(save_path,
                  dpi=400,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None,
                  fmt='png')
    plt.close()
    if os.path.exists(save_path):
        mprint('File {} has created!'.format(save_path))
    else:
        mprint('Error: File {} has no created!'.format(save_path))


def save_graph_lc_2(save_path, JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta, time_interval):
    draw_LC_2(JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta, time_interval)
    plt.savefig(save_path,
                  dpi=400,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None,
                  fmt='png')
    plt.close()
    if os.path.exists(save_path):
        mprint('File {} has created!'.format(save_path))
    else:
        mprint('Error: File {} has no created!'.format(save_path))


def save_aperture(xy, photo, save_path):
    x_star, y_star, x_bckgrnd_in, y_bckgrnd_in, x_bckgrnd_out, y_bckgrnd_out = circle_coord(photo)
    draw_aperture(xy, photo.x, photo.y, x_star, y_star, x_bckgrnd_in, y_bckgrnd_in, x_bckgrnd_out, y_bckgrnd_out)
    plt.savefig(save_path,
                  dpi=400,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None,
                  fmt='png')
    plt.close()
    if os.path.exists(save_path):
        mprint('File {} has created!'.format(save_path))
    else:
        mprint('Error: File {} has no created!'.format(save_path))



def view_lc(time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta):
    draw_LC(time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta)
    plt.show()
    plt.close()


def view_full_frame(xy):
    draw_full_frame(xy)
    plt.show()
    plt.close()


def view_aperture(xy, photo):
    x_star, y_star, x_bckgrnd_in, y_bckgrnd_in, x_bckgrnd_out, y_bckgrnd_out = circle_coord(photo)
    draw_aperture(xy, photo.x, photo.y, x_star, y_star, x_bckgrnd_in, y_bckgrnd_in, x_bckgrnd_out, y_bckgrnd_out)
    plt.show()
    plt.close()


def crop_data(file_path, x_0, y_0, crop_value):
    '''
    Вырезает участок изображения размером 2 * crop_value, с центров в точке x_0, y_0.
    :param file_path: путь к фотонному листу.
    :param x_0: координата х центра области выреза.
    :param y_0: координата y центра области выреза.
    :param crop_value: половина величины размера выреза изображения.
    :return: записывает файл с данными из вырезанной области. Возвращет путь к записанному файлу.
    '''
    save_path = file_path[:-6] + '_crop.ascii'
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line[0] == '#':
                data.append(line)
            else:
                temp_1 = line.strip().split()
                if -crop_value <= int(temp_1[1]) - x_0 < crop_value and -crop_value <= int(temp_1[2]) - y_0 < crop_value:
                    data.append(line)

    with open(save_path, 'w') as f:
        for line in data:
            f.write(line)
    mprint('File {} has created!'.format(save_path))
    del data
    Eat()

    return save_path


def gif_animation():
    pass


def view_light_curve():
    pass


def generate_xy(file_path):
    mprint('Downloading image ...')
    n = 512
    xy = np.zeros((n, n))

    with open(file_path, 'r') as f:
        for index, line in enumerate(f):
            if line[0] != '#':
                temp_1 = line.strip().split()
                xy[int(temp_1[2])][int(temp_1[1])] += 1
    mprint('Image uploaded.')

    return xy


# do = {
#     "view": view_full_frame,
#     "--set_aper": set_aperture,
#     "--set_bkgnd": set_background,
#     "--gif_anim": gif_animation,
#     "--view_lc": view_light_curve
# }


def open_new_data(file_path, photo):
    '''
    Загружает данные из фотонных листов. Определяет попадание
    фотона в соответствующие области: звезда, фон, другое.
    :param file_path: путь к фотонному листу.
    :param photo: параметры фотометрии.
    :return: данные с информацией из фотонных листов.
    '''
    flag = True
    time_start = 0
    data = []

    with open(file_path, 'r') as f:
        for line in f:
            if line[0] != '#':
                temp_1 = line.strip().split()
                if flag:
                    time_start = float(temp_1[0]) 
                    flag = False
                src = (float(temp_1[0]) - time_start, int(temp_1[1]), int(temp_1[2]))
                #if 2500 <= src[0] <= 2600:
                src = RAW_Curve(*src)
                data.append(divide_b_tip_one(src, photo))

    return data


@jit
def show_field_foo(one_time):
    d_field = 512
    temp = np.zeros((d_field, d_field))
    for src in one_time:
        #new_x = src.x - x_0_star + d_field2
        #new_y = src.y - y_0_star + d_field2
        temp[src.y, src.x] += 1

    return temp


def show_field(lst_divide_by_time, x_0_star, y_0_star, x_star, y_star, x_bckgrnd_in, y_bckgrnd_in, x_bckgrnd_out, y_bckgrnd_out, signal, bkgnd, time_interval):
    d_field = 512
    for index, one_time in enumerate(lst_divide_by_time):
        temp = show_field_foo(one_time)
        fig = plt.figure()
        plt.imshow(temp)
        plt.plot(x_0_star, y_0_star, linestyle=' ', marker='.', markersize=3, color='blue')
        try:
            plt.plot(x_star, y_star, linestyle=' ', marker='.', markersize=1, color='red',
                     label='$S/N = {:.1f}$'.format(signal[index] / bkgnd[index]))
        except ZeroDivisionError:
            plt.plot(x_star, y_star, linestyle=' ', marker='.', markersize=1, color='red',
                     label='$S/N = error$')
        plt.plot(x_bckgrnd_in, y_bckgrnd_in, linestyle=' ', marker='.', markersize=1, color='green')
        plt.plot(x_bckgrnd_out, y_bckgrnd_out, linestyle=' ', marker='.', markersize=1, color='green')
        plt.colorbar()
        plt.gca().invert_yaxis()
        save_path = 'C:\\Users\\vasil\\Desktop\\2\\{:05}'.format(index)
        plt.gca().set_ylabel('Y', fontsize=14)
        plt.gca().set_xlabel('X', fontsize=14)
        plt.gca().set_title('{:05} $bin = {}\ s$'.format(index, str(time_interval)))
        #plt.legend()
        plt.savefig(save_path)
        #plt.show()
        plt.close()
        del fig
        del temp
        del one_time
        if not index % 25:
            Eat()


def divide_b_tip_one(src, photo):
    '''
    Определяет попадание фотона в соответствующие области: звезда, фон, другое.
    :param src: данные по одному фотону.
    :param photo: параметры фотометрии.
    :return: данные по одному фотону, с добавлением информации об области попадания этого фотона.
    '''
    if math.sqrt((src.x - photo.x) ** 2 + (src.y - photo.y) ** 2) <= photo.radius_star:
        temp = Curve(*src, in_star='star')
    elif photo.radius_star + photo.gap < math.sqrt((src.x - photo.x) ** 2 + (src.y - photo.y) ** 2) <= photo.radius_star + photo.gap + photo.width_bkgnd:
        temp = Curve(*src, in_star='background')
    else:
        temp = Curve(*src, in_star='other')

    return temp


def get_square(data):
    n = 512
    sq_star = 0
    sq_bkgnd = 0
    temp_data = np.zeros((n, n))
    for src in data:
        if not temp_data[src.y][src.x]:
            temp_data[src.y][src.x] = 1
            if src.in_star == 'star':
                sq_star += 1
            elif src.in_star == 'background':
                sq_bkgnd += 1

    ratio = sq_bkgnd / sq_star

    return sq_star, sq_bkgnd, ratio


def divide_by_time(data, time_interval):
    lst_divide_by_time = []
    temp = data[0].time
    temp_lst = []
    count_quantums = []
    counter_quantums = 0
    for src in data:
        if src.time - temp <= time_interval:
            temp_lst.append(src)
            counter_quantums += 1
        else:
            lst_divide_by_time.append(temp_lst)
            temp = src.time
            temp_lst = [src]
            count_quantums.append(counter_quantums)
            counter_quantums = 0

    return lst_divide_by_time, count_quantums


def convert_png_to_gif():
    dir_path = 'C:\\Users\\vasil\\Desktop\\2'
    gif_dir_path = os.path.join(dir_path, 'gif')

    lst_png = os.scandir(dir_path)
    for src in lst_png:
        if src.is_file():
            img = Image.open(src.path)
            gif_name = src.name[:-4] + '.gif'
            img_gif_path = os.path.join(gif_dir_path, gif_name)
            img.save(img_gif_path)
            os.remove(src.path)


def convert_data_to_signal(lst_divide_by_time, ratio):
    time = []
    signal = []
    err_signal = []
    bkgnd = []
    err_bkgnd = []
    for src in lst_divide_by_time:
        temp_time = np.mean(list(map(lambda x: x.time, src)))
        temp_signal = len(list(filter(lambda x: x.in_star == 'star', src)))
        temp_err_signal = math.sqrt(temp_signal)
        temp_bkgnd = len(list(filter(lambda x: x.in_star == 'background', src))) / ratio
        temp_err_bkgnd = math.sqrt(temp_bkgnd)

        time.append(temp_time)
        signal.append(temp_signal)
        bkgnd.append(temp_bkgnd)
        err_signal.append(temp_err_signal)
        err_bkgnd.append(temp_err_bkgnd)

    delta = list(map(lambda x, y: x - y, signal, bkgnd))
    err_delta = list(map(lambda x, y: x + y, err_signal, err_bkgnd))

    return time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta


def get_JD_from_time(time, my_data_JD):
    JD_0 = JDN.get_JD(my_data_JD.year,
                      my_data_JD.month,
                      my_data_JD.day,
                      my_data_JD.hour,
                      my_data_JD.min,
                      my_data_JD.sec)
    JD = [JDN.get_JD(my_data_JD.year,
                      my_data_JD.month,
                      my_data_JD.day,
                      my_data_JD.hour,
                      my_data_JD.min,
                      float(my_data_JD.sec) + t) for t in time]

    return JD


def start_auto():
    path_json = r'C:\Users\vasil\Desktop\qwe.json'
    data = json.load(open(path_json))
    for ph_sheet in data:
        try:
            data_to_auto(ph_sheet)
            mprint('{:-^40}'.format('end'))
        except Exception as e:
            logging.error('{}: {}'.format(e, ph_sheet))

def mprint(s):
    logging.info(s)
    print(s)


def create_gif_animation(filenames, save_path):

    with imageio.get_writer(save_path, mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
            os.remove(filename)

    if os.path.exists(save_path):
        mprint('File {} has created!'.format(save_path))
    else:
        mprint('Error: File {} has no created!'.format(save_path))


def data_to_auto(ph_sheet):
    mprint(re.findall(r'^(\d{6}_\d{2})', ph_sheet['ascii'])[0])
    tmp_ascii_path = os.path.join(ph_sheet['path'], ph_sheet['ascii'])
    tmp_json_path = os.path.join(ph_sheet['path'], ph_sheet['json'])
       
    settings_photometry = Photometry(**json.load(open(tmp_json_path)))

    mprint('>upload file xy')
    xy = generate_xy(tmp_ascii_path)

    mprint('>save aperture')
    save_ap_path = tmp_ascii_path[:-6] + '_aperture.png'
    save_aperture(xy, settings_photometry, save_ap_path)

    mprint('>upload file lc')
    data_for_lc = open_new_data(tmp_ascii_path, settings_photometry)
    sq_star, sq_bkgnd, ratio = get_square(data_for_lc)

    for time_interval in (0.1, 1):
        mprint('time_interval: {}'.format(time_interval))
        mprint('>get lc')
        pre_filename = tmp_ascii_path[:-6] + '_{}sec'.format(time_interval)
        lst_divide_by_time, count_quantums = divide_by_time(data_for_lc, time_interval)
        Eat()
        time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta = convert_data_to_signal(lst_divide_by_time, ratio)
        
        mprint('>save lc')
        save_graph_path = pre_filename + '_lc.png'
        save_graph_path2 = pre_filename + '_lc2.png'
        save_json_path = pre_filename + '_lc.json'
        save_csv_path = pre_filename + '_lc.csv'

        header_JD, data_JD = open_dict_csv('LC_to_UT.txt')
        my_data_JD = None
        for src in data_JD:
            if src['name'] == re.findall(r'^(\d{6}_\d{2})', ph_sheet['ascii'])[0]:
                my_data_JD = JD_data(**src)
                break
        else:
            mprint('There is no the name in the list.')

        if my_data_JD is not None:
            JD = get_JD_from_time(time, my_data_JD)
            header_csv = ['time', 'JD', 'signal', 'err_signal', 'bkgnd', 'err_bkgnd', 'delta', 'err_delta']
            data = list(zip(time, JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta))
            data_csv = []
            for src in data:
                temp = dict(zip(header_csv, src))
                data_csv.append(temp)
        save_graph_lc_1(save_graph_path, JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta)
        save_graph_lc_2(save_graph_path2, JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta, time_interval)
        save_json_photo(save_json_path, settings_photometry)
        write_dict_csv(save_csv_path, header_csv, data_csv, delimiter=',')

    mprint('>Preparation for obtaining images')
    gif_file_path = tmp_ascii_path[:-6] + '.gif'
    time_interval = 10
    lst_divide_by_time, count_quantums = divide_by_time(data_for_lc, time_interval)
    Eat()
    mprint('>save png --> gif')
    time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta = convert_data_to_signal(lst_divide_by_time, ratio)
    x_star, y_star, x_background_in, y_background_in, x_background_out, y_background_out = circle_coord(settings_photometry)
    show_field(lst_divide_by_time, settings_photometry.x, settings_photometry.y, x_star, y_star, x_background_in, y_background_in, x_background_out, y_background_out, signal, bkgnd, time_interval)
    convert_png_to_gif()

    mprint('>create gif animation')
    filenames = list(map(lambda x: x.path, os.scandir('C:\\Users\\vasil\\Desktop\\2\\gif\\')))
    create_gif_animation(filenames, gif_file_path)



def main():
    #file_path = os.path.abspath(input('Enter the path to file with data:\n'))
    path_json = r'C:\Users\vasil\YandexDisk\WorkPlace\Scripts\Git\Nomus\flare_star\settings.json'
    file_path = r'C:\Users\vasil\Desktop\12\130103_00_crop.ascii'
    settings_photometry = None
    data_for_lc = None
    time = None
    signal = None
    err_signal = None
    bkgnd = None
    err_bkgnd = None
    delta = None
    err_delta = None
    xy = None
    ans = ''
    time_interval = 10
    pre_filename = file_path[:-6] + '_{}sec'.format(time_interval)
    lst_divide_by_time = None
    while True:
        ans = input('Your ans: ').strip()
        if ans == 'upload file xy':
            # file_path = os.path.abspath(input('Enter path to photonic sheet: ').strip())
            xy = generate_xy(file_path)
        elif ans == 'upload file lc':
            if settings_photometry:
                data_for_lc = open_new_data(file_path, settings_photometry)
            else:
                mprint('Load the photometry data (enter: "set aperture").')
            # file_path = os.path.abspath(input('Enter path to photonic sheet: ').strip())
        elif ans == 'get lc':
            if data_for_lc is not None:
                sq_star, sq_bkgnd, ratio = get_square(data_for_lc)
                lst_divide_by_time, count_quantums = divide_by_time(data_for_lc, time_interval)
                del data_for_lc
                Eat()
                time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta = convert_data_to_signal(lst_divide_by_time, ratio)
                view_lc(time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta)
            else:
                mprint('Load the data to build a light curve. (enter: "upload file lc").')
        elif ans == 'view lc':
            if all([time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta]):
                view_lc(time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta)
            else:
                mprint('Load the data to build a light curve. (enter: "get lc").')
        elif ans == 'save lc':
            if all([time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta]):
                save_graph_path = pre_filename + '_lc.png'
                save_graph_path2 = pre_filename + '_lc2.png'
                save_json_path = pre_filename + '_lc.json'
                save_csv_path = pre_filename + '_lc.csv'
                mprint(pre_filename)
                ans = input('Enter the file name to calc Julian date (without name - no): ').strip()
                if ans != 'no':
                    header_JD, data_JD = open_dict_csv('LC_to_UT.txt')
                    my_data_JD = None
                    for src in data_JD:
                        if src['name'] == ans:
                            my_data_JD = JD_data(**src)
                            break
                    else:
                        mprint('There is no the name in the list. (enter: "save lc")')

                    if my_data_JD is not None:
                        JD = get_JD_from_time(time, my_data_JD)
                        header_csv = ['time', 'JD', 'signal', 'err_signal', 'bkgnd', 'err_bkgnd', 'delta', 'err_delta']
                        data = list(zip(time, JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta))
                        data_csv = []
                        for src in data:
                            temp = dict(zip(header_csv, src))
                            data_csv.append(temp)
                    save_graph_lc_1(save_graph_path, JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta)
                    save_graph_lc_2(save_graph_path2, JD, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta, time_interval)
                    save_json_photo(save_json_path, settings_photometry)
                    write_dict_csv(save_csv_path, header_csv, data_csv, delimiter=',')



                else:
                    header_csv = ['time', 'signal', 'err_signal', 'bkgnd', 'err_bkgnd', 'delta', 'err_delta']
                    data = list(zip(time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta))
                    data_csv = []
                    for src in data:
                        temp = dict(zip(header_csv, src))
                        data_csv.append(temp)
                    save_graph_lc(save_graph_path, time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta)
                    save_json_photo(save_json_path, settings_photometry)
                    write_dict_csv(save_csv_path, header_csv, data_csv, delimiter=',')
            else:
                mprint('No data to save. (enter: "get lc").')
        elif ans == 'save aperture':
            if settings_photometry:
                if xy is not None:
                    save_ap_path = pre_filename + '_aperture.png'
                    save_aperture(xy, settings_photometry, save_ap_path)
                else:
                    mprint('No data to save. (enter: "upload file xy").')
            else:
                mprint('No data to save. (enter: "set aperture").')
        elif ans == 'view':
            view_full_frame(xy)
        elif ans == 'save png':
            if settings_photometry is not None:
                if lst_divide_by_time is not None:
                    x_star, y_star, x_background_in, y_background_in, x_background_out, y_background_out = circle_coord(settings_photometry)
                    show_field(lst_divide_by_time, settings_photometry.x, settings_photometry.y, x_star, y_star, x_background_in, y_background_in, x_background_out, y_background_out, signal, bkgnd, time_interval)
                    convert_png_to_gif()
                else:
                    mprint("There is no division of the photon sheet by time (enter: 'upload file lc' -> 'get lc')")
            else:
                mprint('No data of photometry. (enter: "set aperture").')
        elif ans == 'view aperture':
            if not settings_photometry:
                settings_photometry = set_photometry(xy, path_json)
            view_aperture(xy, settings_photometry)
        elif ans == 'view photometry':
            mprint(settings_photometry)
        elif ans == 'save photometry':
            save_json_photo("settings.json", settings_photometry)
        elif ans == 'set aperture':
            settings_photometry = set_photometry(xy, path_json)
        elif ans == 'crop':
            while True:
                crop_ans = input('Take the values of the center of the region to be cut from the file (y/n):').strip()
                if crop_ans == 'q':
                    break
                elif crop_ans == 'y':
                    crop_value = int(input('Enter the crop_value: ').strip())
                    if settings_photometry:
                        crop_file_path = crop_data(file_path, settings_photometry.x, settings_photometry.y, crop_value)
                    else:
                        mprint('Load the photometry data (enter: "set aperture").')
                    break
                elif crop_ans == 'n':
                    crop_value = int(input('Enter the crop_value: ').strip())
                    x_0 = int(input('Enter x0: ').strip())
                    y_0 = int(input('Enter y0: ').strip())
                    crop_file_path = crop_data(file_path, x_0, y_0, crop_value)
                    break
                else:
                    mprint("Wrong key specified")

            
        elif ans == 'q':
            break
        else:
            mprint("Wrong key specified")


if __name__ == '__main__':
    # main()
    start_auto()