#-*- coding: utf-8 -*-

import math, numpy as np
import matplotlib.pyplot as plt
import pdb


"""
Plotting
"""
def decorate(ax, title=None, xlab=None, ylab=None, legend=None, xlim=None, ylim=None, min_yspan=None):
    ax.xaxis.grid(color='k', linestyle='-', linewidth=0.2)
    ax.yaxis.grid(color='k', linestyle='-', linewidth=0.2)
    if xlab: ax.xaxis.set_label_text(xlab)
    if ylab: ax.yaxis.set_label_text(ylab)
    if title: ax.set_title(title, {'fontsize': 20 })
    if legend is not None:
        if legend == True: ax.legend(loc='best')
        else: ax.legend(legend, loc='best')
    if xlim is not None: ax.set_xlim(xlim[0], xlim[1])
    if ylim is not None: ax.set_ylim(ylim[0], ylim[1])
    if min_yspan is not None: ensure_yspan(ax, min_yspan)

def ensure_yspan(ax, yspan):
    ymin, ymax = ax.get_ylim()
    if ymax-ymin < yspan:
        ym =  (ymin+ymax)/2
        ax.set_ylim(ym-yspan/2, ym+yspan/2)

def prepare_fig(fig=None, window_title=None, figsize=(20.48, 10.24), margins=None):
    if fig is None:
        fig = plt.figure(figsize=figsize)
    else:
        plt.figure(fig.number)
    if margins is not None:
        left, bottom, right, top, wspace, hspace = margins
        fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top,
                            hspace=hspace, wspace=wspace)
    #pdb.set_trace()
    if window_title is not None:
         fig.canvas.set_window_title(window_title)
    return fig


def plot_in_grid(time, plots, ncol, figure=None, window_title="None", legend=None, filename=None,
                 margins=None, lw=2.):
    nrow = math.ceil(len(plots)/float(ncol))
    figsize = (10.24*ncol, 2.56*nrow)
    figure = prepare_fig(figure, window_title, figsize=figsize, margins=margins)
    for i, (title, ylab, min_yspan, data) in enumerate(plots):
        ax = figure.add_subplot(nrow, ncol, i+1)
        ax.plot(time, data, linewidth=lw)
        decorate(ax, title=title, ylab=ylab, min_yspan=min_yspan)
    if legend is not None:
        ax.legend(legend, loc='best')
    #save_if(filename)
    return figure
