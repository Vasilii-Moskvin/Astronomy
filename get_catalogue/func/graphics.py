import matplotlib.ticker
import matplotlib.pyplot as plt
import os.path
import numpy as np
from seaborn import heatmap
from func.StarFromCSV import StarFromCSV
from func.work_with_csv import get_column


def save_correlogram(dir_path, rhosns, filt):
    heatmap(rhosns, annot=True, fmt=".3f", linewidths=.5)
    plt.yticks(rotation=0)
    plt.title('Correlogram of stars with the least error ({} filter)\n'.format(filt), fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=16)
    savename = os.path.join(dir_path, 'corr_{}'.format(filt))
    plt.savefig(savename,
                  dpi=None,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None,
                  fmt='svg')
    plt.close()


def save_mnk(x, y, y_mnk, legendtxt, xlab, ylab, tlt, savename):
    plt.plot(x, y_mnk, label=legendtxt)
    plt.legend(loc=2)
    plt.plot(x, y, 'ro')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(tlt)
    plt.grid(True)
    plt.savefig(savename,
                  dpi=None,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None)
    plt.close()


def save_regress(save_path, x, y, y_mnk, tlt, x_lbl, y_lbl, legend_lbl):
    plt.plot(x, y, linestyle=' ', marker='.', color='grey')
    plt.plot(x, y_mnk, color='red', label=legend_lbl)

    legend1 = plt.legend(loc=2, fontsize=16, prop={'size':14}, frameon=True, framealpha=0.5)
    frame = legend1.get_frame()
    frame.set_color('white')
    frame.set_edgecolor('black')

    plt.axis([9, 16, 9, 16])
    formatter = matplotlib.ticker.FormatStrFormatter("%.0f")
    plt.xticks(fontsize=14, fontname='Times')
    plt.yticks(fontsize=14, fontname='Times')
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.gca().set_ylabel(y_lbl, fontsize=14, fontname='Times')
    plt.gca().set_xlabel(x_lbl, fontsize=14, fontname='Times')
    plt.gca().set_title(tlt, fontsize=16, fontname='Times')
    plt.grid(color='grey', linestyle='--', linewidth=0.5)

    plt.savefig(save_path,
                  dpi=None,
                  facecolor='w',
                  edgecolor='w',
                  orientation='portrait',
                  parpertype=None,
                  format=None,
                  transparent=False,
                  bbox_inches=None,
                  pad_inches=0.1,
                  frameon=None,
                  fmt='svg')
    plt.close()


def save_regress_line(dir_path, reducted_catalogue, ABC_by_filt):
    for filt in StarFromCSV.filts:
        x = get_column(reducted_catalogue, '{}_avr_mag'.format(filt))
        y = get_column(reducted_catalogue, '{}_mag_cat'.format(filt))

        A = np.vstack([x, np.ones(len(x))]).T
        k, b = np.linalg.lstsq(A, y)[0]

        y_mnk = [k * i + b for i in x]

        save_name = os.path.join(dir_path, 'regress_{}'.format(filt))
        xlab = 'Instrumental catalogue, mag'
        ylab = 'Catalogue, mag'
        tlt = 'Catalohue vs instrumental magnitudes ({} filter)'.format(filt)
        if ABC_by_filt[filt].b != 0:
            filt2 = StarFromCSV.choose_filter(filt)
            legendtxt = '${}_c={:.3}*{}_i{:+.3}*({}_i-{}_i){:+.3}$\n$y={:.3}*x{:+.3}$'.format(filt, ABC_by_filt[filt].a, filt, ABC_by_filt[filt].b, filt2, filt, ABC_by_filt[filt].c, k, b)
        else:
            legendtxt = '${}_c={:.3}*{}_i{:+.3}$\n$y={:.3}*x{:+.3}$'.format(filt, ABC_by_filt[filt].a, filt, ABC_by_filt[filt].c, k, b)
        save_mnk(x, y, y_mnk, legendtxt, xlab, ylab, tlt, save_name)
