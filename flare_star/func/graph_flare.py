import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import numpy as np


def draw_full_frame(xy):
    sns.heatmap(xy)
    plt.gca().invert_yaxis()
    x_start, x_end = plt.gca().get_xlim()
    y_start, y_end = plt.gca().get_ylim()
    x = np.arange(x_start, x_end + 1, 50)
    y = np.arange(y_start, y_end + 1, 50)
    plt.gca().xaxis.set_ticks(x)
    plt.gca().yaxis.set_ticks(y)
    plt.gca().xaxis.set_major_formatter(ticker.FormatStrFormatter('%0d'))
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%0d'))
    plt.gca().set_ylabel('Y', fontsize=14)
    plt.gca().set_xlabel('X', fontsize=14)


def draw_aperture(xy, x_0_star, y_0_star, x_star, y_star, x_bckgrnd_in, y_bckgrnd_in, x_bckgrnd_out, y_bckgrnd_out):
    draw_full_frame(xy)
    plt.plot(x_star, y_star, linestyle=' ', marker='.', markersize=1, color='red')
    plt.plot(x_bckgrnd_in, y_bckgrnd_in, linestyle=' ', marker='.', markersize=1, color='green')
    plt.plot(x_bckgrnd_out, y_bckgrnd_out, linestyle=' ', marker='.', markersize=1, color='green')
    plt.plot(x_0_star, y_0_star, linestyle=' ', marker='.', markersize=3, color='blue')


def draw_LC(time, signal, err_signal, bkgnd, err_bkgnd, delta, err_delta):#, rec_bkgnd):
    number = 3
    f, axarr = plt.subplots(number, sharex=True)
    #axarr[0].plot(time, signal, linestyle=' ', marker='.', markersize=1, color='red')
    axarr[0].errorbar(time, signal, linestyle=' ', yerr=err_signal, marker='', markersize=1, color='gray')
    axarr[0].plot(time, signal, linestyle=' ', marker='.', markersize=1, color='red')
    #axarr[0].plot(time, bkgnd, linestyle=' ', marker='.', markersize=1, color='green')
    axarr[0].set_title('Light curve of full frame')
    #set_ylim = axarr[0].get_ylim()
    #axarr[0].set_ylim(-20, 130)
    #axarr[0].set_ylabel('count/sec', fontsize=14)

    axarr[1].errorbar(time, delta, linestyle=' ', yerr=err_delta, marker='', markersize=1, color='gray')
    axarr[1].plot(time, delta, linestyle=' ', marker='.', markersize=1, color='blue')
    axarr[1].set_title('Light curve of star')
    axarr[1].set_ylabel('Intensity, count/sec', fontsize=14)
    #axarr[1].set_ylim(set_ylim)

    axarr[2].errorbar(time, bkgnd, linestyle=' ', yerr=err_bkgnd, marker='', markersize=1, color='gray')
    axarr[2].plot(time, bkgnd, linestyle=' ', marker='.', markersize=1, color='green')
    #axarr[2].plot(time, rec_bkgnd, linestyle=' ', marker='.', markersize=1, color='blue')
    axarr[2].set_title('Background')
    #axarr[2].set_ylabel('Star and background, count/sec', fontsize=14)
    axarr[2].set_xlabel('time, s', fontsize=14)
    #axarr[2].set_ylim(set_ylim)

    # axarr[3].plot(time, list(map(lambda x, y: x / y, signal, bkgnd)), linestyle=' ', marker='.', markersize=1, color='green')
    #axarr[3].plot(time, rec_bkgnd, linestyle=' ', marker='.', markersize=1, color='blue')
    # axarr[3].set_title('Rate')
    # axarr[2].set_ylabel('Star and background, count/sec', fontsize=14)
    # axarr[3].set_xlabel('time, s', fontsize=14)
    #axarr[3].set_ylim(set_ylim)

    for index in range(number):
        axarr[index].grid(color='grey', linestyle='--', linewidth=0.5)
        