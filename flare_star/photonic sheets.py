from func.graph_flare import draw_full_frame, draw_aperture
from collections import namedtuple
import math
import numpy as np
import json
import matplotlib.pyplot as plt

Photometry = namedtuple('Photometry', ('x', 'y', 'radius_star', 'width_bkgnd'))


def set_photometry(xy, path_json):
    ans = ''
    while True:
        ans = input(r'Photometry settings (file\auto): ').strip()
        if ans == 'file':
            photometry_settings = Photometry(**json.load(open(path_json)))
            break
        elif ans == 'auto aperture':
            default_settings = Photometry(**json.load(open(path_json)))
            photometry_settings = calc_aperture(xy, default_settings.x, default_settings.y)
            break
        elif ans == 'auto':
            photometry_settings = None
            break
        elif ans == 'q':
            photometry_settings = None
            break
        else:
            print("Wrong key specified. Enter file or auto")

    return photometry_settings


def calc_aperture(xy, x0, y0):
    pass


def circle_coord(photo):
    fi = [i for i in range(360)]
    x_star = [photo.radius_star * math.cos(i) + photo.x for i in fi]
    y_star = [photo.radius_star * math.sin(i) + photo.y for i in fi]
    x_background_inside = [(photo.radius_star + photo.width_bkgnd) * math.cos(i) + photo.x for i in fi]
    y_background_insede = [(photo.radius_star + photo.width_bkgnd) * math.sin(i) + photo.y for i in fi]
    x_background_outside = [(photo.radius_star + photo.width_bkgnd) * math.cos(i) + photo.x for i in fi]
    y_background_outsede = [(photo.radius_star + photo.width_bkgnd) * math.sin(i) + photo.y for i in fi]

    return x_star, y_star, x_background_inside, y_background_insede, x_background_outside, y_background_outsede



def view_full_frame(xy):
    draw_full_frame(xy)
    plt.show()
    plt.close()


def view_aperture(xy, photo):
    x_star, y_star, x_bckgrnd_in, y_bckgrnd_in, x_bckgrnd_out, y_bckgrnd_out = circle_coord(photo)
    draw_aperture(xy, photo.x, photo.y, x_star, y_star, x_bckgrnd_in, y_bckgrnd_in, x_bckgrnd_out, y_bckgrnd_out)
    plt.show()
    plt.close()


def set_aperture():
    pass


def set_background():
    pass


def gif_animation():
    pass


def view_light_curve():
    pass





def generate_xy(file_path):
    print('Downloading image ...')
    n = 512
    xy = np.zeros((n, n))

    with open(file_path, 'r') as f:
        for index, line in enumerate(f):
            if line[0] != '#':
                temp_1 = line.strip().split()
                xy[int(temp_1[2])][int(temp_1[1])] += 1
    print('Image uploaded.')

    return xy


do = {
    "view": view_full_frame,
    "--set_aper": set_aperture,
    "--set_bkgnd": set_background,
    "--gif_anim": gif_animation,
    "--view_lc": view_light_curve
}

 

def main():
    #file_path = os.path.abspath(input('Enter the path to file with data:\n'))
    file_path = r'C:\Users\vasil\Desktop\temp\090128_03_crop.ascii'
    path_json = r'C:\Users\vasil\YandexDisk\WorkPlace\Scripts\Git\Nomus\flare_star\settings.json'
    xy = generate_xy(file_path)
    settings_photometry = set_photometry(xy, path_json)
    ans = ''
    while True:
        ans = input('Your ans: ').strip()
        if ans == 'view':
            view_full_frame(xy)
        elif ans == 'view aperture':
            view_aperture(xy, settings_photometry)
        elif ans == 'q':
            break
        else:
            print("Wrong key specified")
    print(settings_photometry)


if __name__ == '__main__':
    main()