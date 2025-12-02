#-*- coding: utf-8 -*-

import math, numpy as np
import matplotlib.pyplot as plt
import pdb

#import pat3.frames as p3_fr # FIXME... cross include :(
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
         fig.canvas.manager.set_window_title(window_title)
    return fig


def plot_in_grid(time, plots, ncol, figure=None, axes=None, window_title="None", legend=None, filename=None,
                 margins=None, lw=2., extra_rows=0):
    nrow = math.ceil(len(plots)/float(ncol)) + extra_rows
    figsize = (10.24*ncol, 2.56*nrow)
    figure = prepare_fig(figure, window_title, figsize=figsize, margins=margins)
    axes = figure.subplots(nrow, ncol, sharex=True) if axes is None else axes
    for i, (title, ylab, min_yspan, data) in enumerate(plots):
        ax = axes.flatten()[i]
        ax.plot(time, data, linewidth=lw)
        decorate(ax, title=title, ylab=ylab, min_yspan=min_yspan)
    if legend is not None:
        ax.legend(legend, loc='best')
    #save_if(filename)
    return figure, axes



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

    
def plot_3D_traj(ref_traj=None, X=None, fig=None, ax=None, val=None, tc='b', rtc='g', rgtc='k', tlab='aircraft trajectory'):
    fig = fig if fig is not None else plt.figure()
    ax = ax if ax is not None else fig.add_subplot(111, projection='3d')
    if ref_traj is not None:
        pts = np.asarray(ref_traj.get_points())
        xs, ys, zs = pts[:,0], pts[:,1], pts[:,2] 
        ax.plot(xs, ys, zs, color=rtc, label='ref trajectory')
        ax.plot(xs, ys, np.zeros(len(xs)), linestyle='--', color=rgtc, linewidth=0.5, label='ground track')
        # https://stackoverflow.com/questions/36013063/what-is-the-purpose-of-meshgrid-in-python-numpy
        if 1:
            verts = []
            for i in range(1, len(pts)):
                verts.append([(xs[i-1], ys[i-1], 0), (xs[i-1], ys[i-1], zs[i-1]),
                              (xs[i], ys[i], zs[i]), (xs[i], ys[i], 0)])
            mesh = Poly3DCollection(verts, alpha=0.2, linewidth=0, facecolor=[0.5, 0.5, 1])
            ax.add_collection3d(mesh)

    if X is not None:
        if val is None:
            ax.plot(X[:,0], X[:,1], X[:,2], color=tc, label=tlab)
        else:
            #pdb.set_trace()
            _skip=5
            #ac_pos_ned = X[:, p3_fr.SixDOFEuclidianEuler.slice_pos].reshape(-1, 1, 3)
            ac_pos_ned = X[::_skip, :3].reshape(-1, 1, 3)
            segments_ = np.concatenate([ac_pos_ned[:-1], ac_pos_ned[1:]], axis=1)
            #norm = plt.Normalize(_val.min(), _val.max())
            norm = plt.Normalize(0., 1.5)
            lc = mpl_toolkits.mplot3d.art3d.Line3DCollection(segments_, cmap=plt.get_cmap('viridis'), norm=norm)
            lc.set_array(-val[::_skip]) 
            lc.set_linewidth(2)
            ax.add_collection3d(lc)#, zs=z, zdir='z')

    ax.set_xlabel('North')
    ax.set_ylabel('East')
    ax.set_zlabel('Down')        
    set_3D_axes_equal()
    plt.legend(loc='lower left')#'best')

def plot_3D_wind(atm, x0=0, xspan=50, dx=5., h0=-10, hspan=150, dh=10., figure=None, ax=None):
    fig = figure if figure is not None else plt.figure()
    ax = ax if ax is not None else fig.add_subplot(111, projection='3d')
    # Make the grid
    x_r, y_r, z_r = np.arange(-(xspan-x0), xspan-x0, dx), np.arange(-50., 50., 5), np.arange(h0, h0+hspan, dh)
    x, y, z = np.meshgrid(x_r, y_r, z_r)
    wx, wy, wz = np.meshgrid(x_r, y_r, z_r)
    nx, ny, nz = x.shape
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                pos = [x[ix, iy, iz], y[ix, iy, iz], -z[ix, iy, iz]]  # we plot z axis up
                wx[ix, iy, iz], wy[ix, iy, iz], wz[ix, iy, iz] = atm.get_wind(pos, t=0)
    w = np.vstack((wx[np.newaxis],wy[np.newaxis],wz[np.newaxis]))
    w_norm = np.linalg.norm(w, axis = 0)
    c = (w_norm.ravel() - w_norm.min()) / w_norm.ptp()
    c = np.concatenate((c, np.repeat(c, 2)))
    # Colormap
    #c = plt.cm.hsv(c)
    cmap=plt.get_cmap('viridis')
    c = cmap(c)
    #pdb.set_trace()
    #q = ax.quiver(x, y, z, wx, wy, wz)
    #q = ax.quiver(x, y, z, wx, wy, wz, length=1., cmap='Reds
    q = ax.quiver(x, y, z, wx, wy, -wz, colors=c, length=1., lw=3, normalize=True)
    #q = ax.quiver(x, y, z, wx, wy, wz, length=1.1, cmap='Reds', lw=2, normalize=True)
    #q.set_array(np.random.rand(np.prod(x.shape)))
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')



def plot_slice_wind(atm, xmax=50, dx=5., h0=-10, h1=150, dh=2.):
    y=0
    xlist, zlist = np.arange(-xmax, xmax, dx), np.arange(h0, h1, dh)
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
    ax.axis('equal')

#
# vertical wind slice
#
def plot_slice_wind_nu(atm, n0=-50, n1=50, dn=5., e0=0., h0=-10, h1=150, dh=2., t0=0., zdir=-1.,
                       show_quiver=True, contour_wz=False, show_color_bar=True,
                       title=None,
                       figure=None, ax=None):
    fig = figure if figure is not None else plt.figure()
    ax = ax if ax is not None else fig.add_subplot(111)
    _nlist, _ulist = np.arange(n0, n1, dn), np.arange(h0, h1, dh)
    _n, _u = np.meshgrid(_nlist, _ulist, indexing='ij')
    #wx, wz = np.meshgrid(_nlist, _ulist)
    _wind = np.zeros((len(_nlist), len(_ulist), 3))
    for _in in range(len(_nlist)):
        for _iu in range(len(_ulist)):
            pos_ned = [_n[_in, _iu], e0, zdir*_u[_in, _iu]]  # we plot with z axis up
            #wx[_in, _iu], _, wz[_in, _iu] = atm.get_wind_ned(pos_ned, t=0)
            _wind[_in, _iu] = atm.get_wind_ned(pos_ned, t=0)

    if show_quiver:
        #pdb.set_trace()
        q = ax.quiver(_nlist, _ulist, _wind[:,:,0].T, -_wind[:,:,2].T, angles='uv', units='width')
    #if show_quiver: q = ax.quiver(xlist, zlist, wx, zdir*wz, units='width')

    # wind is 'ned', we plot 'en' but contourf wants xy meshgrid indexing (as opposed to ij)
    _d = -_wind[:,:,2].T if contour_wz else np.linalg.norm(_wind, axis=2).T
    cp = ax.contourf(_nlist, _ulist, _d, alpha=0.4)
    
    #cp = ax.contourf(x, z, zdir*wz, alpha=0.4)
    #if use_wx: cp = ax.contourf(x, z, zdir*wx, alpha=0.4) # FIX ME : this can be made more generic for x-y-z
    if show_color_bar:
        cbar = fig.colorbar(cp)
        lab = 'wz in m/s (up)' if contour_wz else 'wind in m/s'
        cbar.ax.set_ylabel(lab, rotation=270)#; cbar.ax.set_xlabel('thermal')
    title = f'NorthUp slice at e:{e0}m t:{t0}s' if title is None else title
    decorate(ax, title=title, xlab='North [$m$]', ylab='Altitude [m] (positive up)', legend=None, xlim=None, ylim=None, min_yspan=None)
    ax.axis('equal')
    return fig, ax


def plot_slice_wind_eu(atm, n0=0, e0=0, e1=100, de=5., h0=-10, h1=150, dh=2., t0=0., zdir=-1.,
                       show_quiver=True, contour_wz=False, show_color_bar=True,
                       title=None,
                       figure=None, ax=None):
    fig = figure if figure is not None else plt.figure()
    ax = ax if ax is not None else fig.add_subplot(111)
    _elist, _ulist = np.arange(e0, e1, de), np.arange(h0, h1, dh)
    _e, _u = np.meshgrid(_elist, _ulist, indexing='ij')
    _wind = np.zeros((len(_elist), len(_ulist), 3))
    for ie in range(len(_elist)):
        for iu in range(len(_ulist)):
            pos_ned = [n0, _e[ie, iu], zdir*_u[ie, iu]]  # we plot with z axis up
            _wind[ie, iu] = atm.get_wind_ned(pos_ned, t=0)
    if show_quiver:
        q = ax.quiver(_elist, _ulist, _wind[:,:,1].T, -_wind[:,:,2].T, angles='uv', units='width')
    # wind is 'ned', we plot 'en' but contourf wants xy meshgrid indexing (as opposed to ij)
    _d = -_wind[:,:,2].T if contour_wz else np.linalg.norm(_wind, axis=2).T
    cp = ax.contourf(_elist, _ulist, _d, alpha=0.4)
    if show_color_bar:
        cbar = fig.colorbar(cp)
        lab = 'wz in m/s (up)' if contour_wz else 'wind in m/s'
        cbar.ax.set_ylabel(lab, rotation=270)#; cbar.ax.set_xlabel('thermal')
    title = f'EastUp slice at n:{n0}m t:{t0}s' if title is None else title
    decorate(ax, title=title, xlab='east in m', ylab='h in m (positive up)', legend=None, xlim=None, ylim=None, min_yspan=None)
    ax.axis('equal')
    return fig, ax


#
# horizontal wind slice
#
def plot_slice_wind_ne(atm, n0=-100, n1=100, dn=5., e0=-100., e1=100, de=5, h0=0., t0=0.,
                       show_quiver=False, contour_wz=False, show_color_bar=False,
                       title=None,
                       figure=None, ax=None):
    fig = figure if figure is not None else plt.figure()
    ax = ax if ax is not None else fig.add_subplot(111)
    _nlist, _elist = np.arange(n0, n1, dn), np.arange(e0, e1, de)
    _n, _e = np.meshgrid(_nlist, _elist, indexing='ij')
    _wind = np.zeros((len(_nlist), len(_elist), 3))
    for _in in range(_wind.shape[0]):
        for _ie in range(_wind.shape[1]):
            pos_ned = [_n[_in, _ie], _e[_in, _ie], -h0]
            _wind[_in, _ie] = atm.get_wind_ned(pos_ned, t=0)
    
    if show_quiver:
        q = ax.quiver(_elist, _nlist, _wind[:,:,1], _wind[:,:,0], angles='uv', units='width')

    # wind is 'ned', we plot 'en' but contourf wants xy meshgrid indexing (as opposed to ij)
    _d = -_wind[:,:,2] if contour_wz else np.linalg.norm(_wind, axis=2)
    cp = ax.contourf(_elist, _nlist, _d, alpha=0.4) 

    if show_color_bar:
        cbar = fig.colorbar(cp)
        lab = 'wz in m/s (up)' if contour_wz else 'wind in m/s'
        cbar.ax.set_ylabel(lab, rotation=270)#; cbar.ax.set_xlabel('wind speed')
    title = f'EastNorth slice at {h0}m {t0}s' if title is None else title
    decorate(ax, title=title, xlab='east in m', ylab='north in m', legend=None, xlim=None, ylim=None, min_yspan=None)
    ax.axis('equal')
    return figure, ax
