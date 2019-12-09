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
                 margins=None, lw=2., extra_rows=0):
    nrow = math.ceil(len(plots)/float(ncol)) + extra_rows
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



'''
3D
'''
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def set_3D_axes_equal(ax=None):
    '''
    Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''
    if ax is None: ax = plt.gca()

    x_limits, y_limits, z_limits = ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()

    x_range, x_middle = abs(x_limits[1] - x_limits[0]), np.mean(x_limits)
    y_range, y_middle = abs(y_limits[1] - y_limits[0]), np.mean(y_limits)
    z_range, z_middle = abs(z_limits[1] - z_limits[0]), np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    
def plot_3D_traj(ref_traj=None, X=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if ref_traj is not None:
        pts = ref_traj.get_points()
        xs, ys, zs = pts[:,0], pts[:,1], pts[:,2] 
        ax.plot(xs, ys, zs, color='g', label='ref trajectory')
        ax.plot(xs, ys, np.zeros(len(xs)), linestyle='--', color='k', linewidth=0.5, label='ground track')
        # https://stackoverflow.com/questions/36013063/what-is-the-purpose-of-meshgrid-in-python-numpy
        if 1:
            verts = []
            for i in range(1, len(pts)):
                verts.append([(xs[i-1], ys[i-1], 0), (xs[i-1], ys[i-1], zs[i-1]),
                              (xs[i], ys[i], zs[i]), (xs[i], ys[i], 0)])
            mesh = Poly3DCollection(verts, alpha=0.2, linewidth=0, facecolor=[0.5, 0.5, 1])
            ax.add_collection3d(mesh)

    if X is not None:
        ax.plot(X[:,0], X[:,1], X[:,2], color='b', label='aircraft trajectory')

    ax.set_xlabel('North')
    ax.set_ylabel('East')
    ax.set_zlabel('Down')        
    set_3D_axes_equal()
    plt.legend(loc='best')

def plot_3D_wind(atm):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # Make the grid
    x_r, y_r, z_r = np.arange(-50., 50., 5), np.arange(-50., 50., 5), np.arange(0., 100, 10)
    x, y, z = np.meshgrid(x_r, y_r, z_r)
    wx, wy, wz = np.meshgrid(x_r, y_r, z_r)
    nx, ny, nz = x.shape
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                pos = [x[ix, iy, iz], y[ix, iy, iz], z[ix, iy, iz]]
                wx[ix, iy, iz], wy[ix, iy, iz], wz[ix, iy, iz] = atm.get_wind(pos, t=0)
    w = np.vstack((wx[np.newaxis],wy[np.newaxis],wz[np.newaxis]))
    w_norm = np.linalg.norm(w, axis = 0)
    c = (w_norm.ravel() - w_norm.min()) / w_norm.ptp()
    c = np.concatenate((c, np.repeat(c, 2)))
    # Colormap
    c = plt.cm.hsv(c)
    #pdb.set_trace()
    #q = ax.quiver(x, y, z, wx, wy, wz)
    #q = ax.quiver(x, y, z, wx, wy, wz, length=1., cmap='Reds
    q = ax.quiver(x, y, z, wx, wy, wz, colors=c, length=1., lw=3, normalize=True)
    #q = ax.quiver(x, y, z, wx, wy, wz, length=1.1, cmap='Reds', lw=2, normalize=True)
    #q.set_array(np.random.rand(np.prod(x.shape)))
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

def plot_slice_wind(atm, xmax=50, dx=5.):
    y=0
    xlist, zlist = np.arange(-xmax, xmax, dx), np.arange(-10., 150, 2)
    x, z = np.meshgrid(xlist, zlist)
    wz = np.zeros_like(z)
    nx, nz = x.shape
    for ix in range(x.shape[0]):
        for iz in range(x.shape[1]):
            wz[ix, iz] = -atm.get_wind([x[ix, iz], y, -z[ix, iz]], t=0)[2]
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(x, z, wz)
    fig.colorbar(cp)
    ax.set_title('Filled Contours Plot')
