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