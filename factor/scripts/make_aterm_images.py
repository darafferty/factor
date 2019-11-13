#! /usr/bin/env python
"""
Script to make a-term images from solutions
"""
from losoto.h5parm import h5parm
from losoto.operations.directionscreen import _calc_piercepoint
from losoto.operations.plotscreen import _calculate_screen
import lsmtool
from numpy import kron, concatenate, newaxis
from numpy.linalg import pinv, norm
import os
import numpy as np
from factor.lib import miscellaneous as misc
from astropy.io import fits as pyfits
from astropy import wcs
from shapely.geometry import Point, Polygon
from scipy.spatial import Voronoi
import shapely.geometry
import shapely.ops
import scipy.ndimage as ndimage
from lofarpipe.support.data_map import DataMap, DataProduct


def save_frame(screen, fitted_phase1, residuals,  x, y, k, sindx,
    root_dir, prestr, midRA, midDec, order, outQueue):
    """
    Saves screen images as FITS files

    Parameters
    ----------
    screen : array
        Image of screen values
    fitted_phase1 : array
        Array of fitted phase values
    residuals : array
        Array of phase residuals at the piercepoints
    weights : array
        Array of weights at the piercepoints
    x : array
        Array of piercepoint x locations
    y : array
        Array of piercepoint y locations
    k : int
        Time index
    order : int
        order of screen
    is_phase : bool
        True if screens are phase screens

    """
    if not 'matplotlib' in sys.modules:
        import matplotlib as mpl
        mpl.rc('figure.subplot',left=0.05, bottom=0.05, right=0.95, top=0.95,wspace=0.22, hspace=0.22 )
        mpl.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt # after setting "Agg" to speed up
    from losoto.operations.stationscreen import _makeWCS, _circ_chi2
    import numpy as np
    try:
        try:
            from astropy.visualization.wcsaxes import WCSAxes
            hasWCSaxes = True
        except:
            from wcsaxes import WCSAxes
            hasWCSaxes = True
    except:
        hasWCSaxes = False
    from matplotlib.colors import LinearSegmentedColormap


    fig = plt.figure(figsize=(6,6))

    # Set colormap
    if is_phase:
        cmap = _phase_cm()
    else:
        cmap = plt.cm.jet
    sm = plt.cm.ScalarMappable(cmap=cmap,
        norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []

    if is_image_plane and hasWCSaxes:
        wcs = _makeWCS(midRA, midDec)
        ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
    else:
        plt.gca().set_aspect('equal')
        ax = plt.gca()

    s = []
    c = []
    xf = []
    yf = []
    weights = np.array(weights, dtype=float)
    nonflagged = np.where(weights > 0.0)
    for j in range(fitted_phase1.shape[0]):
        if weights[j] > 0.0:
            s.append(max(20, 200*np.sqrt(weights[j]/np.median(weights[nonflagged]))))
        else:
            s.append(120)
            xf.append(x[j])
            yf.append(y[j])
        c.append(sm.to_rgba(fitted_phase1[j]))

    if is_image_plane:
        min_x = np.min(x)
        max_x = np.max(x)
        min_y = np.min(y)
        max_y = np.max(y)
        extent_x = max_x - min_x
        extent_y = max_y - min_y
        lower = [min_x - 0.1 * extent_x, min_y - 0.1 * extent_y]
        upper = [max_x + 0.1 * extent_x, max_y + 0.1 * extent_y]
        Nx = screen.shape[0]
        pix_per_m = Nx / (upper[0] - lower[0])
        m_per_pix = 1.0 / pix_per_m
        xr = np.arange(lower[0], upper[0], m_per_pix)
        yr = np.arange(lower[1], upper[1], m_per_pix)
        lower = np.array([xr[0], yr[0]])
        upper = np.array([xr[-1], yr[-1]])
    else:
        # convert from m to km
        lower /= 1000.0
        upper /= 1000.0

    im = ax.imshow(screen.transpose([1, 0])[:, :],
        cmap = cmap,
        origin = 'lower',
        interpolation = 'nearest',
        extent = (lower[0], upper[0], lower[1], upper[1]),
        vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(im)
    cbar.set_label('Value', rotation=270)

    ax.scatter(np.array(x), np.array(y), s=np.array(s), c=np.array(c), alpha=0.7, cmap=cmap, vmin=vmin, vmax=vmax, edgecolor='black')
    if len(xf) > 0:
        ax.scatter(xf, yf, s=120, c='k', marker='x')
    if show_source_names:
        labels = source_names
        for label, xl, yl in zip(labels, x, y):
            plt.annotate(
                label,
                xy = (xl, yl), xytext = (-2, 2),
                textcoords = 'offset points', ha = 'right', va = 'bottom')

    nsrcs = np.where(weights > 0.0)[0].size
    if is_phase:
        redchi2 =  _circ_chi2(residuals, weights) / (nsrcs-order)
    else:
        redchi2 =  np.sum(np.square(residuals) * weights) / (nsrcs-order)
    if sindx >= 0:
        plt.title('Station {0}, Time {1} (red. chi2 = {2:0.3f})'.format(station_names[sindx], k, redchi2))
    else:
        plt.title('Time {0}'.format(k))
    if is_image_plane:
        ax.set_xlim(lower[0], upper[0])
        ax.set_ylim(lower[1], upper[1])
        ax.set_aspect('equal')
        if hasWCSaxes:
            RAAxis = ax.coords['ra']
            RAAxis.set_axislabel('RA', minpad=0.75)
            RAAxis.set_major_formatter('hh:mm:ss')
            DecAxis = ax.coords['dec']
            DecAxis.set_axislabel('Dec', minpad=0.75)
            DecAxis.set_major_formatter('dd:mm:ss')
            ax.coords.grid(color='black', alpha=0.5, linestyle='solid')
            plt.xlabel("RA")
            plt.ylabel("Dec")
        else:
            plt.xlabel("RA (arb. units)")
            plt.ylabel("Dec (arb. units)")
    else:
        # Reverse the axis so that RA coord increases to left
        plt.xlim(upper[0], lower[0])
        plt.ylim(lower[1], upper[1])
        plt.xlabel('Projected Distance East-West (km)')
        plt.ylabel('Projected Distance North-South (km)')
    if sindx >= 0:
        plt.savefig(root_dir + '/' + prestr + '_station%0.4i' % sindx + '_frame%0.4i.png' % k, bbox_inches='tight')
    else:
        plt.savefig(root_dir + '/' + prestr + '_frame%0.4i.png' % k, bbox_inches='tight')
    plt.close(fig)


def make_screen_images(pp, inscreen, inresiduals, weights, station_names,
    station_positions, source_names, times, height, station_order, beta_val,
    r_0, prefix='frame_', remove_gradient=True, show_source_names=False,
    min_val=None, max_val=None, is_phase=False, midRA=0.0, midDec=0.0, ncpu=0):
    """
    Makes a-term images of screens

    Parameters
    ----------
    pp : array
        Array of piercepoint locations
    inscreen : array
        Array of screen values at the piercepoints with order [dir, time, ant]
    residuals : array
        Array of screen residuals at the piercepoints with order [dir, time, ant]
    weights : array
        Array of weights for each piercepoint with order [dir, time, ant]
    source_names: array
        Array of source names
    times : array
        Array of times
    height : float
        Height of screen (m)
    station_order : list
        List of order of screens per station (e.g., number of KL base vectors to keep)
    r_0 : float
        Scale size of phase fluctuations
    beta_val : float
        Power-law index for phase structure function
    prefix : str
        Prefix for output file names
    remove_gradient : bool
        Fit and remove a gradient from each screen
    show_source_names : bool
        Label sources on screen plots
    min_val : float
        Minimum value for plot range
    max_val : float
        Maximum value for plot range
    is_phase : bool
        Input screen is a phase screen
    midRA : float
        RA for WCS reference in degrees
    midDec : float
        Dec for WCS reference in degrees
    ncpu : int
        Number of CPUs to use

    """
    # input check
    root_dir = os.path.dirname(prefix)
    if root_dir == '':
        root_dir = './'
    prestr = os.path.basename(prefix) + 'screen'
    try:
        os.makedirs(root_dir)
    except OSError:
        pass

    N_stations = 1 # screens are single-station screens
    N_sources = len(source_names)
    N_times = len(times)
    N_piercepoints = N_sources * N_stations
    xp, yp, zp = station_positions[0, :] # use first station
    east = np.array([-yp, xp, 0])
    east = east / norm(east)

    north = np.array([-xp, -yp, (xp*xp + yp*yp)/zp])
    north = north / norm(north)

    up = np.array([xp, yp, zp])
    up = up / norm(up)

    T = concatenate([east[:, newaxis], north[:, newaxis]], axis=1)

    is_image_plane = True # pierce points are image plane coords
    pp1_0 = pp[:, 0:2]
    pp1_1 = pp[:, 0:2]

    max_xy = np.amax(pp1_0, axis=0) - np.amin(pp1_0, axis=0)
    max_xy_1 = np.amax(pp1_1, axis=0) - np.amin(pp1_1, axis=0)
    if max_xy_1[0] > max_xy[0]:
        max_xy[0] = max_xy_1[0]
    if max_xy_1[1] > max_xy[1]:
        max_xy[1] = max_xy_1[1]

    min_xy = np.array([0.0, 0.0])
    extent = max_xy - min_xy
    lower = min_xy - 0.1 * extent
    upper = max_xy + 0.1 * extent
    im_extent_m = upper - lower
    fitted_phase1 = np.zeros((N_piercepoints, N_times))

    Nx = 60 # set approximate number of pixels in screen
    pix_per_m = Nx / im_extent_m[0]
    m_per_pix = 1.0 / pix_per_m
    xr = np.arange(lower[0], upper[0], m_per_pix)
    yr = np.arange(lower[1], upper[1], m_per_pix)
    Nx = len(xr)
    Ny = len(yr)
    lower = np.array([xr[0], yr[0]])
    upper = np.array([xr[-1], yr[-1]])

    x = np.zeros((N_times, N_piercepoints)) # plot x pos of piercepoints
    y = np.zeros((N_times, N_piercepoints)) # plot y pos of piercepoints
    screen = np.zeros((Nx, Ny, N_times))

    for sindx in range(station_positions.shape[0]):
        logging.info('Calculating screen images...')
        residuals = inresiduals[:, :, sindx, newaxis].transpose([0, 2, 1]).reshape(N_piercepoints, N_times)
        mpm = multiprocManager(ncpu, _calculate_screen)
        for k in range(N_times):
            mpm.put([inscreen[:, k, sindx], residuals[:, k], pp,
                N_piercepoints, k, east, north, up, T, Nx, Ny, sindx, height,
                beta_val, r_0, is_phase])
        mpm.wait()
        for (k, ft, scr, xa, ya) in mpm.get():
            screen[:, :, k] = scr
            fitted_phase1[:, k] = ft
            if is_image_plane:
                x[k, :] = xa
                y[k, :] = ya
            else:
                x[k, :] = xa - np.amin(xa) # remove offsets for each time slot
                y[k, :] = ya - np.amin(ya)

        # Normalize piercepoint locations to extent calculated above
        if not is_image_plane:
            x *= extent[0]
            y *= extent[1]

        if min_val is None:
            vmin = np.min([np.amin(screen), np.amin(fitted_phase1)])
        else:
            vmin = min_val
        if max_val is None:
            vmax = np.max([np.amax(screen), np.amax(fitted_phase1)])
        else:
            vmax = max_val

    logging.info('Save screen images...')
    for sindx in range(station_positions.shape[0]):
        for k in range(N_times):
            mpm.put([screen[:, :, k], fitted_phase1[:, k], residuals[:, k],
            weights[:, k, sindx], x[k, :], y[k, :], k, lower, upper, vmin, vmax,
            source_names, show_source_names, station_names, sindx, root_dir,
            prestr, is_image_plane, midRA, midDec, station_order[0, k, sindx], is_phase])


def gaussian_fcn(g, x1, x2):
    """Evaluate Gaussian on the given grid.

    Parameters
    ----------
    x1, x2: grid (as produced by numpy.mgrid f.e.)
    g: Gaussian object or list of Gaussian paramters
    """
    from math import radians, sin, cos

    A, C1, C2, S1, S2, Th = g
    fwsig=2.35482
    S1 = S1/fwsig
    S2 = S2/fwsig
    Th = Th + 90.0 # Define theta = 0 on x-axis

    th = radians(Th)
    cs = cos(th)
    sn = sin(th)

    f1 = ((x1-C1)*cs + (x2-C2)*sn)/S1
    f2 = (-(x1-C1)*sn + (x2-C2)*cs)/S2

    return A*np.exp(-(f1*f1 + f2*f2)/2)


def guassian_image(A, x, y, xsize, ysize, gsize):
    """Makes an image of a Gaussian

    Parameters
    ----------
    A : peak value
    x, y : pixel coords of Gaussian center
    xsize, ysize : size of image in pixels
    gsize : size as FWHM in pixels
    """
    im = np.zeros((ysize, xsize))
    g = [A, y, x, gsize, gsize, 0.0]
    bbox = np.s_[0:ysize, 0:xsize]
    x_ax, y_ax = np.mgrid[bbox]
    gimg = gaussian_fcn(g, x_ax, y_ax)
    ind = np.where(np.abs(gimg) > abs(A)/3.0)
    if len(ind[0]) > 0:
        im[ind[0], ind[1]] += gimg[ind]

    return im


def main(h5parmfile, soltabname, outroot, bounds_deg, bounds_mid_deg, skymodel,
         solsetname='sol000', ressoltabname='', padding_fraction=1.4, cellsize_deg=0.1,
         smooth_deg=0, gsize_deg=0, mapfile_dir='.', filename='aterm_images.mapfile',
         time_avg_factor=1, fasth5parm=None):
    """
    Make a-term FITS images

    Parameters
    ----------
    h5parmfile : str
        Filename of h5parm
    soltabname : str
        Soltab containing solutions or screen fit
    outroot : str
        Root of filename of output FITS file (root+'_0.fits')
    bounds_deg : list
        List of [maxRA, minDec, minRA, maxDec] for image bounds
    bounds_mid_deg : list
        List of [RA, Dec] for midpoint of image bounds
    skymodel : str
        Filename of calibration sky model (needed for patch positions)
    solsetname : str, optional
        Name of solset
    ressoltabname : str, optional
        Soltab containing the screen residuals
    padding_fraction : float, optional
        Fraction of total size to pad with (e.g., 0.2 => 20% padding all around)
    cellsize_deg : float, optional
        Cellsize of output image
    smooth_deg : float, optional
        Size of smoothing kernel in degrees to apply
    gsize_deg : float, optional
        FWHM in degrees of Gaussian to add at patch locations (to enforce that
        solutions at these locations are exactly equal to those in the h5parm)
    mapfile_dir : str, optional
        Directory in which to store output mapfile
    filename : str, optional
        Filename of output mapfile
    time_avg_factor : int, optional
        Averaging factor in time for fast-phase corrections
    fasth5parm : str, optional
        Filename of fast-phase h5parm to be added together with input h5parm

    Returns
    -------
    result : dict
        Dict with list of FITS files
    """
    # Read in solutions
    H = h5parm(h5parmfile)
    solset = H.getSolset(solsetname)
    if 'gain' in soltabname:
        soltab = solset.getSoltab(soltabname.replace('gain', 'amplitude'))
        soltab_ph = solset.getSoltab(soltabname.replace('gain', 'phase'))
    else:
        soltab = solset.getSoltab(soltabname)

    if type(bounds_deg) is str:
        bounds_deg = [float(f.strip()) for f in bounds_deg.strip('[]').split(';')]
    if type(bounds_mid_deg) is str:
        bounds_mid_deg = [float(f.strip()) for f in bounds_mid_deg.strip('[]').split(';')]
    if padding_fraction is not None:
        padding_fraction = float(padding_fraction)
        padding_ra = (bounds_deg[2] - bounds_deg[0]) * (padding_fraction - 1.0)
        padding_dec = (bounds_deg[3] - bounds_deg[1]) * (padding_fraction - 1.0)
        bounds_deg[0] -= padding_ra
        bounds_deg[1] -= padding_dec
        bounds_deg[2] += padding_ra
        bounds_deg[3] += padding_dec
    cellsize_deg = float(cellsize_deg)
    gsize_deg = float(gsize_deg)
    gsize_pix = gsize_deg / cellsize_deg
    smooth_deg = float(smooth_deg)
    smooth_pix = smooth_deg / cellsize_deg
    time_avg_factor = int(time_avg_factor)

    if 'screen' in soltab.getType():
        # Handle screen solutions
        screen_type = soltab.getType()
        logging.info('Using input {0} soltab {1}'.format(screen_type, soltab.name))

        # Get values from soltabs
        solset = soltab.getSolset()
        if resSoltab is '':
            try:
                # Look for residual soltab assuming standard naming conventions
                ressoltab = solset.getSoltab(soltab.name+'resid')
            except:
                logging.error('Could not find the soltab with associated screen residuals. '
                    'Please specify it with the "resSoltab" argument.')
                return 1
        else:
            ressoltab = solset.getSoltab(resSoltab)
        logging.info('Using input screen residual soltab: {}'.format(ressoltab.name))
        screen = np.array(soltab.val)
        weights = np.array(soltab.weight)
        residuals = np.array(ressoltab.val)
        times = np.array(soltab.time)
        orders = np.array(ressoltab.weight)
        axis_names = soltab.getAxesNames()
        if len(screen.shape) > 3:
            # remove degenerate freq and pol axes
            if 'freq' in axis_names:
                freq_ind = axis_names.index('freq')
                screen = np.squeeze(screen, axis=freq_ind)
                weights = np.squeeze(weights, axis=freq_ind)
                residuals = np.squeeze(residuals, axis=freq_ind)
                orders = np.squeeze(orders, axis=freq_ind)
                axis_names.remove('freq')
            if 'pol' in axis_names:
                pol_ind = axis_names.index('pol')
                screen = np.squeeze(screen, axis=pol_ind)
                weights = np.squeeze(weights, axis=pol_ind)
                residuals = np.squeeze(residuals, axis=pol_ind)
                orders = np.squeeze(orders, axis=pol_ind)
                axis_names.remove('pol')

        # Rearrange to get order [dir, time, ant]
        dir_ind = axis_names.index('dir')
        time_ind = axis_names.index('time')
        ant_ind = axis_names.index('ant')
        screen = screen.transpose([dir_ind, time_ind, ant_ind])
        weights = weights.transpose([dir_ind, time_ind, ant_ind])
        residuals = residuals.transpose([dir_ind, time_ind, ant_ind])
        orders = orders.transpose([dir_ind, time_ind, ant_ind])

        # Collect station and source names and positions and times, making sure
        # that they are ordered correctly.
        source_names = soltab.dir[:]
        source_dict = solset.getSou()
        source_positions = []
        for source in source_names:
            source_positions.append(source_dict[source])
        station_names = soltab.ant
        station_dict = solset.getAnt()
        station_positions = []
        for station in station_names:
            station_positions.append(station_dict[station])
        height = soltab.obj._v_attrs['height']
        beta_val = soltab.obj._v_attrs['beta']
        r_0 = soltab.obj._v_attrs['r_0']
        pp = soltab.obj.piercepoint[:]
        if height == 0.0:
            midRA = soltab.obj._v_attrs['midra']
            midDec = soltab.obj._v_attrs['middec']
        else:
            midRA = 0.0
            midDec = 0.0

        if (minZ == 0 and maxZ == 0):
            min_val = None
            max_val = None
        else:
            min_val = minZ
            max_val = maxZ

        _make_screen_plots(pp, screen, residuals, weights, np.array(station_names),
            np.array(station_positions), np.array(source_names), times,
            height, orders, beta_val, r_0, prefix=prefix,
            remove_gradient=remove_gradient, show_source_names=show_source_names, min_val=min_val,
            max_val=max_val, is_phase=is_phase, midRA=midRA, midDec=midDec, ncpu=ncpu)

    else:
        # Handle non-screen solutions using Voronoi tessellation + smoothing
        if 'amplitude' in soltab.getType():
            # complexgain
            vals = soltab.val
            vals_ph = soltab_ph.val
        else:
            # scalarphase -> set amplitudes to unity
            vals_ph = soltab.val
            vals = np.ones_like(vals_ph)
        times = soltab.time
        freqs = soltab.freq
        ants = soltab.ant
        axis_names = soltab.getAxesNames()
        source_names = soltab.dir[:]

        # Load fast-phase solutions if needed and combine with slow gains by interpolating
        # the slow gains to the fast time grid
        if 'amplitude' in soltab.getType() and fasth5parm is not None:
            H_fast = h5parm(fasth5parm)
            solset_fast = H_fast.getSolset('sol000')
            soltab_fast = solset_fast.getSoltab('phase000')
            times_fast = soltab_fast.time

            # Interpolate the slow gains to the fast times
            axis_names = soltab.getAxesNames()  # assume both have same axes
            time_ind = axis_names.index('time')
            f = si.interp1d(times, vals, axis=time_ind, kind='nearest', fill_value='extrapolate')
            vals = f(times_fast)
            f = si.interp1d(times, vals_ph, axis=time_ind, kind='nearest', fill_value='extrapolate')
            vals_ph = f(times_fast) + soltab_fast.val
            times = times_fast

        # Make blank output FITS file (type does not matter at this point)
        midRA = bounds_mid_deg[0]
        midDec = bounds_mid_deg[1]
        temp_image = outroot + '.tmp'
        imsize = (bounds_deg[3] - bounds_deg[1])  # deg
        imsize = int(imsize / cellsize_deg)  # pix
        misc.make_template_image(temp_image, midRA, midDec, ximsize=imsize,
                                 yimsize=imsize, cellsize_deg=cellsize_deg, freqs=freqs,
                                 times=[0.0], antennas=soltab.ant, aterm_type='tec')
        hdu = pyfits.open(temp_image, memmap=False)
        data = hdu[0].data
        w = wcs.WCS(hdu[0].header)
        RAind = w.axis_type_names.index('RA')
        Decind = w.axis_type_names.index('DEC')

        # Get x, y coords for directions in pixels. We use the input calibration sky
        # model for this, as the patch positions written to the h5parm file by DPPP may
        # be different
        skymod = lsmtool.load(skymodel)
        source_dict = skymod.getPatchPositions()
        source_positions = []
        for source in source_names:
            radecpos = source_dict[source.strip('[]')]
            source_positions.append([radecpos[0].value, radecpos[1].value])
        source_positions = np.array(source_positions)
        ra_deg = source_positions.T[0]
        dec_deg = source_positions.T[1]
        xy = []
        for RAvert, Decvert in zip(ra_deg, dec_deg):
            ra_dec = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])
            ra_dec[0][RAind] = RAvert
            ra_dec[0][Decind] = Decvert
            xy.append((w.wcs_world2pix(ra_dec, 0)[0][RAind], w.wcs_world2pix(ra_dec, 0)[0][Decind]))

        # Get boundary for imaging region in pixels
        ra_dec = np.array([[0.0, 0.0, 0.0, 0.0, 0.0]])
        ra_dec[0][RAind] = bounds_deg[0]
        ra_dec[0][Decind] = bounds_deg[1]
        field_minxy = (w.wcs_world2pix(ra_dec, 0)[0][RAind], w.wcs_world2pix(ra_dec, 0)[0][Decind])
        ra_dec[0][RAind] = bounds_deg[2]
        ra_dec[0][Decind] = bounds_deg[3]
        field_maxxy = (w.wcs_world2pix(ra_dec, 0)[0][RAind], w.wcs_world2pix(ra_dec, 0)[0][Decind])
        field_poly = Polygon([[field_minxy[0], field_maxxy[0]],
                             [field_minxy[0], field_maxxy[1]],
                             [field_maxxy[0], field_maxxy[1]],
                             [field_maxxy[0], field_minxy[1]]])

        # Generate array of outer points used to constrain the facets
        nouter = 64
        means = np.ones((nouter, 2)) * np.array(xy).mean(axis=0)
        offsets = []
        angles = [np.pi/(nouter/2.0)*i for i in range(0, nouter)]
        for ang in angles:
            offsets.append([np.cos(ang), np.sin(ang)])
        radius = 2.0*np.sqrt( (field_maxxy[0]-field_minxy[0])**2 + (field_maxxy[1]-field_minxy[1])**2 )
        scale_offsets = radius * np.array(offsets)
        outer_box = means + scale_offsets

        # Tessellate and clip
        points_all = np.vstack([xy, outer_box])
        vor = Voronoi(points_all)
        lines = [
            shapely.geometry.LineString(vor.vertices[line])
            for line in vor.ridge_vertices
            if -1 not in line
        ]
        polygons = [poly for poly in shapely.ops.polygonize(lines)]

        # Index polygons to directions
        ind = []
        for i, xypos in enumerate(xy):
            for poly in polygons:
                if poly.contains(Point(xypos)):
                    poly.index = i

        # Rasterize the polygons to an array, with the value being equal to the
        # polygon's index+1
        data_template = np.ones(data[0, 0, 0, :, :].shape)
        data_rasertize_template = np.zeros(data[0, 0, 0, :, :].shape)
        for poly in polygons:
            verts_xy = poly.exterior.xy
            verts = []
            for x, y in zip(verts_xy[0], verts_xy[1]):
                verts.append((x, y))
            poly_raster = misc.rasterize(verts, data_template.copy()) * (poly.index+1)
            filled = np.where(poly_raster > 0)
            data_rasertize_template[filled] = poly_raster[filled]

        # Identify any gaps in time (frequency gaps are not allowed), as we need to
        # output a separate FITS file for each time chunk
        delta_times = times[1:] - times[:-1]  # time at center of solution interval
        timewidth = np.min(delta_times)
        gaps = np.where(delta_times > timewidth*1.2)
        gaps_ind = gaps[0] + 1
        gaps_ind = np.append(gaps_ind, np.array([len(times)]))

        # Add additional breaks to gaps_ind to keep memory use within that available
        # From experience, making a (30, 46, 62, 4, 146, 146) aterm image needs around
        # 30 GB of memory
        if soltab.getType() == 'tec':
            max_ntimes = 15 * 46 * 4
        else:
            max_ntimes = 15
        # TODO: adjust max_ntimes depending on available memory and time_avg_factor
        check_gaps = True
        while check_gaps:
            check_gaps = False
            g_start = 0
            gaps_ind_copy = gaps_ind.copy()
            for gnum, g_stop in enumerate(gaps_ind_copy):
                if g_stop - g_start > max_ntimes:
                    new_gap = g_start + int((g_stop - g_start) / 2)
                    gaps_ind = np.insert(gaps_ind, gnum, np.array([new_gap]))
                    check_gaps = True
                    break
                g_start = g_stop

        if soltab.getType() == 'tec':
            # TEC solutions
            # input data are [time, ant, dir, freq]
            # output data are [RA, DEC, ANTENNA, FREQ, TIME].T
            # Now loop over stations, frequencies, and times and fill in the correct
            # values
            outfiles = []
            g_start = 0
            for gnum, g_stop in enumerate(gaps_ind):
                outfile = '{0}_{1}.fits'.format(outroot, gnum)
                misc.make_template_image(temp_image, midRA, midDec, ximsize=imsize,
                                         yimsize=imsize, cellsize_deg=cellsize_deg,
                                         times=times[g_start:g_stop],
                                         freqs=freqs, antennas=soltab.ant,
                                         aterm_type='tec')
                hdu = pyfits.open(temp_image, memmap=False)
                data = hdu[0].data
                w = wcs.WCS(hdu[0].header)
                for t, time in enumerate(times[g_start:g_stop]):
                    for f, freq in enumerate(freqs):
                        for s, stat in enumerate(ants):
                            for poly in polygons:
                                ind = np.where(data_rasertize_template == poly.index+1)
                                data[t, f, s, ind[0], ind[1]] = vals[t+g_start, s, poly.index, 0]

                            # Smooth if desired
                            if smooth_pix > 0:
                                data[t, f, s, :, :] = ndimage.gaussian_filter(data[t, f, s, :, :],
                                                      sigma=(smooth_pix, smooth_pix), order=0)

                            # Add Gaussians at patch positions if desired
                            if gsize_pix > 0:
                                for i, (x, y) in enumerate(xy):
                                    # Only do this if patch is inside the region of interest
                                    if int(x) >= 0 and int(x) < data.shape[4] and int(y) >= 0 and int(y) < data.shape[3]:
                                        A = vals[t+g_start, s, i, 0] - data[t, f, s, int(y), int(x)]
                                        data[t, f, s, :, :] += guassian_image(A, x, y, data.shape[4],
                                                               data.shape[3], gsize_pix)
                g_start = g_stop

                # Write FITS file
                hdu[0].data = data
                hdu.writeto(outfile, overwrite=True)
                outfiles.append(outfile)
                os.remove(temp_image)

            map_out = DataMap([])
            map_out.data.append(DataProduct('localhost', ','.join(outfiles), False))
            fileid = os.path.join(mapfile_dir, filename)
            map_out.save(fileid)

            return {'aterm_images': ','.join(outfiles)}
        else:
            # Gain solutions
            # input data are [time, freq, ant, dir, pol] for slow gains (complexgain)
            # and [time, freq, ant, dir] for fast (non-tec) phases (scalarphase)
            # output data are [RA, DEC, MATRIX, ANTENNA, FREQ, TIME].T
            # Now loop over stations, frequencies, and times and fill in the correct
            # matrix values (matrix dimension has 4 elements: real XX, imaginary XX,
            # real YY and imaginary YY)
            outfiles = []
            g_start = 0
            for gnum, g_stop in enumerate(gaps_ind):
                outfile = '{0}_{1}.fits'.format(outroot, gnum)
                misc.make_template_image(temp_image, midRA, midDec, ximsize=imsize,
                                         yimsize=imsize, cellsize_deg=cellsize_deg,
                                         times=times[g_start:g_stop],
                                         freqs=freqs, antennas=soltab.ant,
                                         aterm_type='gain')
                hdu = pyfits.open(temp_image, memmap=False)
                data = hdu[0].data
                w = wcs.WCS(hdu[0].header)
                for t, time in enumerate(times[g_start:g_stop]):
                    for f, freq in enumerate(freqs):
                        for s, stat in enumerate(ants):
                            for p, poly in enumerate(polygons):
                                ind = np.where(data_rasertize_template == poly.index+1)
                                if 'pol' in axis_names:
                                    val_amp_xx = vals[t+g_start, f, s, poly.index, 0]
                                    val_amp_yy = vals[t+g_start, f, s, poly.index, 1]
                                    val_phase_xx = vals_ph[t+g_start, f, s, poly.index, 0]
                                    val_phase_yy = vals_ph[t+g_start, f, s, poly.index, 1]
                                else:
                                    val_amp_xx = vals[t+g_start, f, s, poly.index]
                                    val_amp_yy = vals[t+g_start, f, s, poly.index]
                                    val_phase_xx = vals_ph[t+g_start, f, s, poly.index]
                                    val_phase_yy = vals_ph[t+g_start, f, s, poly.index]
                                data[t, f, s, 0, ind[0], ind[1]] = val_amp_xx * np.cos(val_phase_xx)
                                data[t, f, s, 2, ind[0], ind[1]] = val_amp_yy * np.cos(val_phase_yy)
                                data[t, f, s, 1, ind[0], ind[1]] = val_amp_xx * np.sin(val_phase_xx)
                                data[t, f, s, 3, ind[0], ind[1]] = val_amp_yy * np.sin(val_phase_yy)

                            # Smooth if desired
                            if smooth_pix > 0:
                                data[t, f, s, :, :, :] = ndimage.gaussian_filter(data[t, f, s, :, :, :], sigma=(0, smooth_pix, smooth_pix), order=0)

                            # Add Gaussians at patch positions if desired
                            if gsize_pix > 0:
                                for i, (x, y) in enumerate(xy):
                                    # Only do this if patch is inside the region of interest
                                    if int(x) >= 0 and int(x) < data.shape[4] and int(y) >= 0 and int(y) < data.shape[3]:
                                        if 'pol' in axis_names:
                                            val_amp_xx = vals[t+g_start, f, s, i, 0]
                                            val_amp_yy = vals[t+g_start, f, s, i, 1]
                                            val_phase_xx = vals_ph[t+g_start, f, s, i, 0]
                                            val_phase_yy = vals_ph[t+g_start, f, s, i, 1]
                                        else:
                                            val_amp_xx = vals[t+g_start, f, s, i]
                                            val_amp_yy = vals[t+g_start, f, s, i]
                                            val_phase_xx = vals_ph[t+g_start, f, s, i]
                                            val_phase_yy = vals_ph[t+g_start, f, s, i]
                                        A = val_amp_xx * np.cos(val_phase_xx) - data[t, f, s, 0, int(y), int(x)]
                                        data[t, f, s, 0, :, :] += guassian_image(A, x, y, data.shape[5], data.shape[4], gsize_pix)
                                        A = val_amp_yy * np.cos(val_phase_yy) - data[t, f, s, 2, int(y), int(x)]
                                        data[t, f, s, 2, :, :] += guassian_image(A, x, y, data.shape[5], data.shape[4], gsize_pix)
                                        A = val_amp_xx * np.sin(val_phase_xx) - data[t, f, s, 1, int(y), int(x)]
                                        data[t, f, s, 1, :, :] += guassian_image(A, x, y, data.shape[5], data.shape[4], gsize_pix)
                                        A = val_amp_yy * np.sin(val_phase_yy) - data[t, f, s, 3, int(y), int(x)]
                                        data[t, f, s, 3, :, :] += guassian_image(A, x, y, data.shape[5], data.shape[4], gsize_pix)
                g_start = g_stop

                # If averaging in time, make a new template image with
                # fewer times and write to that instead
                if time_avg_factor > 1:
                    times_avg = times[g_start:g_stop:time_avg_factor]
                    misc.make_template_image(temp_image+'.avg', midRA, midDec, ximsize=imsize,
                                             yimsize=imsize, cellsize_deg=cellsize_deg,
                                             times=times_avg,
                                             freqs=freqs, antennas=soltab.ant,
                                             aterm_type='gain')
                    hdu = pyfits.open(temp_image+'.avg', memmap=False)
                    data_avg = hdu[0].data

                    # Average
                    for t, time in enumerate(times_avg):
                        incr = min(time_avg_factor, len(times[g_start:g_stop])-t*time_avg_factor)
                        data_avg[t, :, :, :, :, :] = np.nanmean(data[t:t+incr, :, :, :, :, :], axis=0)
                    data = data_avg

                # Ensure there are no NaNs in the images, as WSClean will produced uncorrected,
                # uncleaned images if so
                for t, time in enumerate(times[g_start:g_stop]):
                    for p in range(4):
                        if p % 2:
                            # Imaginary elements
                            nanval = 0.0
                        else:
                            # Real elements
                            nanval = 1.0
                        data[t, :, :, p, :, :][np.isnan(data[t, :, :, p, :, :])] = nanval

                # Write FITS file
                hdu[0].data = data
                hdu.writeto(outfile, overwrite=True)
                outfiles.append(outfile)
                os.remove(temp_image)
                hdu = None
                data = None

            map_out = DataMap([])
            map_out.data.append(DataProduct('localhost', ','.join(outfiles), False))
            fileid = os.path.join(mapfile_dir, filename)
            map_out.save(fileid)

            return {'aterm_images': ','.join(outfiles)}
