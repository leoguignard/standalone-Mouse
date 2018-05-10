#!python
# This file is subject to the terms and conditions defined in
# file 'LICENSE', which is part of this source code package.
# Author: Leo Guignard (guignardl...@AT@...janelia.hhmi.org)

import numpy as np
import os
from IO import imread, imsave, SpatialImage
from multiprocessing import Pool, Manager, cpu_count
import psutil, shutil
from skimage.feature import register_translation
from pyklb import readheader

def get_position(t, X):
    ''' Given a sparse corresponding values,
        computes the linear interpolation between these values.
        Args:
            t: int, queried value
            X: {int: float, }, dictionary of corresponding values
        Returns:
            int, the interpolated value.
    '''
    times = np.array(sorted(X.keys()))
    pos = np.min(np.where(t < times)[0])
    t_b = times[pos-1]
    t_e = times[pos]
    return np.int(np.round(np.linspace(X[t_b], X[t_e], t_e - t_b + 1)[t - t_b]))


def read_mask_and_save(path, t, X, ors, vs = (2., 2., 10.)):
    ''' Reads, resamples, masks and saves an image from klb to inr format
        Args:
            path: string, path of the input klb image
            t: int, time point to read
            X: int, value above which the image data is considered
            ors: 3x1 array_like, original aspect ratio
            vs: 3x1 array_like, new aspect ratio
        Return:
            im: string, path to the resampled masked image
    '''
    out = ('.TMP/' + os.path.basename(path)%(t)).replace('.klb', '.inr')
    if not os.path.exists(out):
        t_for_string = (t,)*p.count('%')
        im = imread(path%t_for_string)
        im.voxelsize = vs

        mask = (slice(get_position(t, X), -1), slice(0, -1), slice(0, -1))
        im[mask] = 0

        if (ors != vs).any():
            im.voxelsize = ors
            imsave(out, im)
            os.system(path_to_bin + 'applyTrsf ' + out +
                               ' ' + out +
                               ' -resize -vs %f %f %f '%tuple(vs))
        else:
            imsave(out, im)

    return out

def read_mask(path, t, X, ors, vs = (2., 2., 10.)):
    ''' Reads, masks and saves an image from klb to inr format
        Args:
            path: string, path of the input klb image
            t: int, time point to read
            X: int, value above which the image data is considered
            ors: 3x1 array_like, original aspect ratio
            vs: 3x1 array_like, new aspect ratio
        Return:
            im: NxMxL array_like, masked image
    '''
    out = ('.TMP/' + os.path.basename(path)%t).replace('.klb', '.inr')
    if not os.path.exists(out):
        t_for_string = (t,)*p.count('%')
        im = imread(path%t_for_string)
        im.voxelsize = vs
        mask = (slice(get_position(t, X), -1), slice(0, -1), slice(0, -1))
        im[mask] = 0
        imsave(out, im)
    else:
        im = imread(out)
    return im

def read_trsf(path):
    ''' Read a transformation from a text file
        Args:
            path: string, path to a transformation
    '''
    f = open(path)
    if f.read()[0] == '(':
        f.close()
        f = open(path)
        lines = f.readlines()[2:-1]
        f.close()
        return -np.array([[float(v) for v in l.split()]  for l in lines])
    else:
        f.close()
        return -np.loadtxt(path)

def produce_trsf(params):
    ''' Given an output path, an image path format, a reference time, two time points,
        two locks on their respective images, an original voxel size and a voxel size,
        compute the transformation that registers together two consecutive in time images.
        The registration is done by cross correlation of the MIP along the 3 main axes.
        This function is meant to be call by multiprocess.Pool
    '''
    p_out, p, r, t1, t2, l1, l2, ors, vs = params
    if not os.path.exists(p_out + 't%03d-%03d.txt'%(t1, t2)) or os.path.exists(p_out + 't%03d-%03d.txt'%(t2, t1)):
        path_format = p
        reference_time = r
        # for the reading and copying of the images, a lock is acquired to avoid conflictual reading/writing
        l1.acquire()
        out1 = read_mask_and_save(path_format, t1, X, ors, vs)
        l1.release()
        l2.acquire()
        out2 = read_mask_and_save(path_format, t2, X, ors, vs)
        l2.release()
        if t1 < reference_time:
            os.system(path_to_bin + 'blockmatching -ref ' + out2 + ' -flo ' + out1 +\
                      ' -trsf-type translation -py-hl 6 -py-ll 2' + \
                      ' -res-voxel-trsf ' + p_out + 't%03d-%03d.txt'%(t1, t2))
        else:
            os.system(path_to_bin + 'blockmatching -ref ' + out1 + ' -flo ' + out2 +\
                      ' -trsf-type translation -py-hl 6 -py-ll 2' + \
                      ' -res-voxel-trsf ' + p_out + 't%03d-%03d.txt'%(t2, t1))

        # os.system(path_to_bin + 'rm ' + out1 + ' ' + out2)

def produce_trsf_simple(params):
    ''' Given an output path, an image path format, a reference time, two time points,
        two locks on their respective images, an original voxel size and a voxel size,
        compute the transformation that registers together two consecutive in time images.
        The registration is done by cross correlation of the MIP along the 3 main axes.
        This function is meant to be call by multiprocess.Pool
    '''
    p_out, p, r, t1, t2, l1, l2, ors, vs = params
    path_format = p
    reference_time = r

    # for the reading and copying of the images, a lock is acquired to avoid conflictual reading/writing
    l1.acquire()
    im1 = read_mask(path_format, t1, X, ors, vs)
    l1.release()
    l2.acquire()
    im2 = read_mask(path_format, t2, X, ors, vs)
    l2.release()
    if t1 < reference_time:
        im1, im2 = im2, im1
        out_mat = p_out + 't%03d-%03d.txt'%(t1, t2)
    else:
        out_mat = p_out + 't%03d-%03d.txt'%(t2, t1)

    # projection of the two images and then 
    # computation of the translation that minimizes
    # the cross correlation of the two images
    final_shift = np.zeros(3)
    im1_proj = np.max(im1, axis = 2)
    im2_proj = np.max(im2, axis = 2)
    common_shape = np.max([im1_proj.shape, im2_proj.shape], axis = 0)
    c_s_image1 = np.zeros(common_shape)
    c_s_image2 = np.zeros(common_shape)
    c_s_image1[:im1_proj.shape[0], :im1_proj.shape[1]] = im1_proj
    c_s_image2[:im2_proj.shape[0], :im2_proj.shape[1]] = im2_proj
    shift1, error, diffphase = register_translation(c_s_image1, c_s_image2)
    final_shift[:-1] += shift1

    im1_proj = np.max(im1, axis = 0)
    im2_proj = np.max(im2, axis = 0)
    common_shape = np.max([im1_proj.shape, im2_proj.shape], axis = 0)
    c_s_image1 = np.zeros(common_shape)
    c_s_image2 = np.zeros(common_shape)
    c_s_image1[:im1_proj.shape[0], :im1_proj.shape[1]] = im1_proj
    c_s_image2[:im2_proj.shape[0], :im2_proj.shape[1]] = im2_proj
    shift2, error, diffphase = register_translation(c_s_image1, c_s_image2)
    final_shift[1:] += shift2

    im1_proj = np.max(im1, axis = 1)
    im2_proj = np.max(im2, axis = 1)
    common_shape = np.max([im1_proj.shape, im2_proj.shape], axis = 0)
    c_s_image1 = np.zeros(common_shape)
    c_s_image2 = np.zeros(common_shape)
    c_s_image1[:im1_proj.shape[0], :im1_proj.shape[1]] = im1_proj
    c_s_image2[:im2_proj.shape[0], :im2_proj.shape[1]] = im2_proj
    shift3, error, diffphase = register_translation(c_s_image1, c_s_image2)
    final_shift[::2] += shift3

    final_shift /= 2

    mat_trsf = np.zeros((4, 4))
    mat_trsf[:-1, -1] = -final_shift
    mat_trsf[np.diag_indices(4)] = 1

    np.savetxt(out_mat, mat_trsf)

def run_produce_trsf(p, r, X, nb_times, trsf_p, tp_list, ors, first_TP = 0, vs = (3., 3., 5.), nb_cpu = 1):
    ''' Parallel processing of the transformations from t to t-1/t-1 to t (depending on t<r)
        The transformation is computed using blockmatching algorithm
        Args:
            p: string, path pattern to the images to register
            r: int, reference time point
            nb_times (not used): int, number of time points on which to apply the transformation
            trsf_p: string, path to the transformation
            tp_list: [int, ], list of time points on which to apply the transformation
            ors: float, original aspect ratio
            first_TP (not used): int, first time point on which to apply the transformation
            vs: float, aspect ratio
            nb_cpy: int, number of cpus to use
    '''
    from time import time
    if not os.path.exists(trsf_p):
        os.makedirs(trsf_p)
    # images has to be copied, since the process is parallelized,
    # locks are necessary to avoid conflictual readings/writings
    m = Manager()
    locks = {l: m.Lock() for l in tp_list}
    mapping = [(trsf_p, p, r, t1, t2, locks[t1], locks[t2], ors, vs)
                for t1, t2 in zip(tp_list[:-1], tp_list[1:])]# range(first_TP, nb_times + first_TP)]
    pool = Pool(processes = nb_cpu)
    tic = time()
    tmp = pool.map(produce_trsf, mapping)
    tac = time()
    pool.close()
    pool.terminate()
    whole_time = tac - tic
    secs = whole_time%60
    whole_time = whole_time//60
    mins = whole_time%60
    hours = whole_time//60
    print '%dh:%dmin:%ds'%(hours, mins, secs)

def run_produce_trsf_simple(p, r, X, nb_times, trsf_p, tp_list, ors, first_TP = 0, vs = (3., 3., 5.), nb_cpu = 1):
    ''' Parallel processing of the transformations from t to t-1/t-1 to t (depending on t<r)
        The transformation is computed using cross correlation of the MIP
        Args:
            p: string, path pattern to the images to register
            r: int, reference time point
            nb_times (not used): int, number of time points on which to apply the transformation
            trsf_p: string, path to the transformation
            tp_list: [int, ], list of time points on which to apply the transformation
            ors: float, original aspect ratio
            first_TP (not used): int, first time point on which to apply the transformation
            vs: float, aspect ratio
            nb_cpy: int, number of cpus to use
    '''
    from time import time
    if not os.path.exists(trsf_p):
        os.mkdir(trsf_p)
    m = Manager()
    locks = {l: m.Lock() for l in tp_list}
    mapping = [(trsf_p, p, r, t1, t2, locks[t1], locks[t2], ors, vs)
                for t1, t2 in zip(tp_list[:-1], tp_list[1:])]
    pool = Pool(processes = nb_cpu)
    tic = time()
    tmp = pool.map(produce_trsf_simple, mapping)
    tac = time()
    pool.close()
    pool.terminate()
    whole_time = tac - tic
    secs = whole_time%60
    whole_time = whole_time//60
    mins = whole_time%60
    hours = whole_time//60
    print '%dh:%dmin:%ds'%(hours, mins, secs)

def compose_trsf(flo_t, ref_t, trsf_p, tp_list):
    ''' Recusrively build the transformation that allows
        to register time `flo_t` onto the frame of time `ref_t`
        assuming that it exists the necessary intermediary transformations
        Args:
            flo_t: int, time of the floating image
            ref_t: int, time of the reference image
            trsf_p: string, path to folder containing the transformations
            tp_list: [int, ], list of time points that have been processed
        Returns:
            out_trsf: string, path to the result composed transformation
    '''
    out_trsf = trsf_p + 't%03d-%03d.txt'%(flo_t, ref_t)
    if not os.path.exists(out_trsf):
        flo_int = tp_list[tp_list.index(flo_t) + np.sign(ref_t - flo_t)]
        # the call is recursive, to build `T_{flo\leftarrow ref}`
        # we need `T_{flo+1\leftarrow ref}` and `T_{flo\leftarrow ref-1}`
        trsf_1 = compose_trsf(flo_int, ref_t, trsf_p, tp_list)
        trsf_2 = compose_trsf(flo_t, flo_int, trsf_p, tp_list)
        os.system(path_to_bin + 'composeTrsf ' + out_trsf + ' -trsfs ' + trsf_2 + ' ' + trsf_1)
    return out_trsf

def get_trsf(folder_name, nb_times, r, tp_list, vs, ors, first_TP = 0, bm = True):
    ''' Given a path to transformations, a list of time points, a refenrence time,
        a voxel size, an original voxel size reads the transformations.
        Args:
            folder_name: string, path the folder containing the transformations, 
                         the trsf files names should be in the form t%03d-%03d.txt
            nb_times (not used): int, number of time points on which to apply the transformation
            r: int, reference time point
            tp_list: [int, ], list of time points on which to apply the transformation
            vs: float, aspect ratio
            ors: float, original aspect ratio
            first_TP (not used): int, first time point on which to apply the transformation
            bm (not used): bool, whether the transformation were computed with blockmatching or not            
        Returns:
            Mx3 array_like: list of all the transformations ordered as in tp_list
    '''
    trsfs = []
    for t in tp_list:#range(first_TP, nb_times + 1):
        if t != r:
            trsf = folder_name + 't%03d-%03d.txt'%(t, r)
            trsfs.append(read_trsf(trsf)[:-1, -1]*(vs/ors))
        else:
            trsfs.append(np.array([0., 0., 0.]))
    return np.array(trsfs)

def print_trsf(folder_name, folder_out, nb_times, r, tp_list, vs, ors, first_TP = 0, bm = True):
    ''' Given a path to transformations, a list of time points, a refenrence time,
        a voxel size, an original voxel size and whether blockmatching were used or not,
        produces a csv files with a translation per line.
        A line in the csv correspond to a transformation at a given time in that form:
        time, d x, d y, d z
        Args:
            folder_name: string, path the folder containing the transformations, 
                         the trsf files names should be in the form t%03d-%03d.txt
            folder_out: folder where the csv file will be written
            nb_times (not used): int, number of time points on which to apply the transformation
            r: int, reference time point
            tp_list: [int, ], list of time points on which to apply the transformation
            vs: float, aspect ratio
            ors: float, original aspect ratio
            first_TP (not used): int, first time point on which to apply the transformation
            bm (not used): bool, whether the transformation were computed with blockmatching or not
    '''
    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
    trsfs = get_trsf(folder_name, nb_times, r, tp_list, vs, ors, first_TP, bm)
    # fig = plt.figure(figsize = (10, 8))
    # ax = fig.add_subplot(111)

    # ax.plot(trsfs[:, 0], label = 'X correction', lw = 3, alpha = .5)
    # ax.plot(trsfs[:, 1], label = 'Y correction', lw = 3, alpha = .5)
    # ax.plot(trsfs[:, 2], label = 'Z correction', lw = 3, alpha = .5)
    # ax.set_ylabel('Displacement [voxel]', fontsize = 25)
    # ax.set_xlabel('time [TP]', fontsize = 25)
    # ax.legend(loc = 'upper left', fontsize = 20, ncol=2)
    # ax.tick_params(labelsize = 20)
    # fig.tight_layout()
    # plt.savefig(folder_out + 'drift_vs_time.pdf')
    # plt.close()
    f = open(folder_out + 'trsf.csv', 'w')
    for t, (x, y, z) in zip(tp_list, trsfs):
        f.write('%d, %f, %f, %f\n'%(t, x, y, z))
    f.close()

def build_and_shift_proj(im, shift, axis):
    ''' Provided an image, and transformation and an axis.
        Performs the shift according to the transformation and the MIP along that axis.
        Args:
            im: array_like, the image to shift
            shift: 3x1 numpy.array, the translation
            axis: int, axis along which to project
        Returns:
            im_out: array_like, the shifted projection along axis
    '''
    im_proj = np.max(im, axis = axis)
    im_out = np.zeros_like(im_proj)
    x_len, y_len = im_out.shape
    shift = np.round(shift)
    if int(shift[0]) >= 0 and int(shift[1]) >= 0:
        im_out[int(shift[0]):, int(shift[1]):] = im_proj[:int(x_len - shift[0]), :int(y_len - shift[1])]
    elif int(shift[0]) >= 0:
        im_out[int(shift[0]):, :int(shift[1])] = im_proj[:int(x_len - shift[0]), int(- shift[1]):]
    elif int(shift[1]) >= 0:
        im_out[:int(shift[0]), int(shift[1]):] = im_proj[int(- shift[0]):, :int(y_len - shift[1])]
    else:
        im_out[:int(shift[0]), :int(shift[1])] = im_proj[int(- shift[0]):, int(- shift[1]):]
    return im_out

def run_shift_proj(params):
    ''' Given the path to an image, the path to a transformation, at time and a reference time,
        build the 'xy', 'xz', 'yz' MIP of the shifted image.
        Meant to be run with multiprocess.Pool
        Returns:
            - list of the three projections and the time of the image
    '''
    p_im, p_trsf, t, r, vs = params
    im = imread(p_im%t)
    if not os.path.exists(p_trsf%(t, r)):
        final_shift = np.zeros((4, 4))
        final_shift[np.diag_indices(4)] = 1
    else:
        final_shift = read_trsf(p_trsf%(t, r))
    im_proj_1 = build_and_shift_proj(im, [final_shift[0, 3], final_shift[1, 3]], 2)
    im_proj_2 = build_and_shift_proj(im, [final_shift[1, 3], final_shift[2, 3]], 0)
    im_proj_3 = build_and_shift_proj(im, [final_shift[0, 3], final_shift[2, 3]], 1)
    return [im_proj_1, im_proj_2, im_proj_3, t]

def get_offset(trsf, init_global_shape, init_shape):
    ''' Gets the resulting offset of an image knowing an initial global shape
        and the initial shape of the image
        Args:
            trsf: 3 x 1 matrix transformation to apply
            init_global_shape: 3 x 1 matrix initial proposed global shape
            init_shape: 3 x 1 matrix shape of the image to transform
        Returns:
            3 x 1 matrix, minimum potential offset
            3 x 1 matrix, maximum potential offset
    '''
    margin = init_global_shape.astype(int) - init_shape
    new_shape = np.zeros((2, 3))
    new_shape[1, :] = np.copy(init_shape)
    new_shape += trsf
    return (new_shape[0,:], new_shape[1,:])


def assemble_projection(registered_folder, tp_list, folder_name, proj = 'xy'):
    ''' Build a "movie" of the MIP so they are easy to read
        Args:
            registered_folder: string, path to the folder that contain the registered images
            tp_list: [int, ], list of time points to process
            folder_name: string, output folder name
            proj: string ('xy', 'xz', 'yz'), projection to build; default 'xy'
    '''
    t_for_string2 = (tp_list[0],)*registered_folder.count('%')
    registered_folder_t = registered_folder%t_for_string2
    pos_klb = registered_folder_t.find('.klb')
    dims = tuple(readheader(registered_folder_t[:pos_klb] + '_' + proj + 'Projection' + 
                          registered_folder_t[pos_klb:])['imagesize_tczyx'][:2:-1])
    out = np.zeros(dims + (len(tp_list),), dtype = np.uint16)
    for i, t in enumerate(tp_list):
        t_for_string2 = (t,)*registered_folder.count('%')
        registered_folder_t = registered_folder%t_for_string2
        pos_klb = registered_folder_t.find('.klb')
        out[:,:,i] = imread(registered_folder_t[:pos_klb] + '_' + proj + 'Projection' + registered_folder_t[pos_klb:])[:,:,0]

    imsave(folder_name + '_out/' + proj + '_proj.klb', SpatialImage(out))

def apply_single_trsf(params):
    ''' Apply a transformation to an image and build the xy, xz, yz MIP.
        Meant to be called in parallel by a multiprocess.Pool
    '''
    im_p, init_trsf, init_global_shape, trsf, final_shape, registered_folder_t, t = params
    curr_state = 'Not even started'
    # Since the process is error prone, we allow to catch preciselly the errors with the try/catch
    try:
        curr_state = 'reading input image.'
        im = imread(im_p)
        # build the final image for time t
        new_image = np.zeros(final_shape[::-1], dtype = im.dtype)
        starting_point = (trsf - init_trsf).astype(int)
        ending_point = starting_point + im.shape
        curr_state = 'applying trsf.'
        new_image[[slice(starting_point[i], ending_point[i]) for i in range(2, -1, -1)]] = im.transpose(2, 1, 0)
        if not os.path.exists(os.path.dirname(registered_folder_t)):
            os.makedirs(os.path.dirname(registered_folder_t))
        curr_state = 'saving image.'
        imsave(registered_folder_t, new_image)

        # build the final MIP for time t
        curr_state = 'building max projections.'
        pos_klb = registered_folder_t.find('.klb')
        projXY = np.max(new_image, axis=0)
        projXY = projXY.reshape(*((1,) + projXY.shape))
        projXZ = np.max(new_image, axis=1)
        projXZ = projXZ.reshape(*((1,) + projXZ.shape))
        projYZ = np.max(new_image, axis=2)
        projYZ = projYZ.reshape(*((1,) + projYZ.shape))

        curr_state = 'saving xy projection.'
        imsave(registered_folder_t[:pos_klb] + '_xyProjection' + registered_folder_t[pos_klb:], projXY)
        curr_state = 'saving xz projection.'
        imsave(registered_folder_t[:pos_klb] + '_xzProjection' + registered_folder_t[pos_klb:], projXZ)
        curr_state = 'saving yz projection.'
        imsave(registered_folder_t[:pos_klb] + '_yzProjection' + registered_folder_t[pos_klb:], projYZ)
    except Exception as e:
        print 'Error for time point %d.'%t
        print 'Happened while ' + curr_state
        print 'Error raised:'
        print e

def applyTrsf_from_run(trsf_p, nb_times, r, tp_list, vs, ors, first_TP, registered_folder, nb_cpu, TP_to_keep):
    ''' Apply transformations from paths to transforamtions
        Args: 
            trsf_p: string, path the folder containing the transformations, 
                    the trsf files names should be in the form t%03d-%03d.txt
            nb_times (not used): int, number of time points on which to apply the transformation
            r: int, reference time point
            tp_list: [int, ], list of time points on which to apply the transformation
            vs: float, aspect ratio
            ors: float, original aspect ratio
            first_TP (not used): int, first time point on which to apply the transformation
            registered_folder: string, folder that will contain the registered images
            nb_cpu: int. number of cpus to use
            TP_to_keep: [int, ], list of time points to treat (if only as subset is necessary)
    '''
    # reads all the transformations
    trsfs = np.round(get_trsf(trsf_p, nb_times, r, tp_list, vs, ors, first_TP, True)).astype(int)
    shapes = []
    trsfs_dict = dict(zip(tp_list, trsfs))

    # reads the shape of each time point in order to compute the minimum bounding box
    for t in tp_list:
        t_for_string = (t,)*p.count('%')
        shapes += [list(readheader(p%t_for_string)['imagesize_tczyx'][:1:-1])]

    # compute the new transformation to apply the right cropping of the image
    init_global_shape = np.max(shapes, axis = 0)
    min_shape, max_shape = [], []
    for i, t in enumerate(trsfs):
        min_s, max_s = get_offset(t, init_global_shape, shapes[i])
        min_shape += [min_s]
        max_shape += [max_s]
    min_shape = np.min(min_shape, axis = 0)
    max_shape = np.max(max_shape, axis = 0)
    final_shape = np.ceil(max_shape - min_shape).astype(int)

    # apply the transformations in parallel
    mapping = []
    if set(TP_to_keep).intersection(tp_list) != set():
        TP_to_do = set(TP_to_keep).intersection(tp_list)
    else:
        TP_to_do = tp_list
    for i, t in enumerate(TP_to_do):
        t_for_string = (t,)*p.count('%')
        t_for_string2 = (t,)*registered_folder.count('%')
        mapping += [[p%t_for_string, min_shape, init_global_shape, trsfs_dict[t], final_shape, registered_folder%t_for_string2, t]]
    pool = Pool(processes = nb_cpu)
    out = pool.map(apply_single_trsf, mapping)
    pool.close()
    pool.terminate()

def read_param_file():
    ''' Reads a csv parameter file.
        Asserts if all the parameters have been introduced
    '''
    p_param = raw_input('Please enter the path to the parameter file/folder:\n')
    p_param = p_param.replace('"', '')
    p_param = p_param.replace("'", '')
    p_param = p_param.replace(" ", '')
    if p_param[-4:] == '.csv':
        f_names = [p_param]
    else:
        f_names = [os.path.join(p_param, f) for f in os.listdir(p_param) if 'param_reg' in f and '.csv' in f and not '~' in f]
    paths = []
    Xs = []
    reference_times = []
    time_list = []
    VS = []
    OR = []
    suffixes = []
    times = []
    registration_methods = []
    registered_folder = []
    comp_trsfs = []
    app_trsfs = []
    TP_to_keep = []
    paths_to_bin = []
    for file_name in f_names:
        f = open(file_name)
        lines = f.readlines()
        f.close()
        param_dict = {}
        i = 0
        nb_lines = len(lines)
        while i < nb_lines:
            l = lines[i]
            split_line = l.split(',')
            param_name = split_line[0]
            if (param_name != 'crop' and  param_name != 'voxel_downsampled' and param_name != 'TP_to_keep' and
                param_name != 'original_voxel_ratio' and param_name != '' and param_name != 'time_to_remove'):
                val = split_line[1]
                if val.isdigit():
                    val = int(val)
                param_dict[param_name] = val
                i += 1
            elif param_name == 'crop':
                Xs_tmp = {}
                while param_name == 'crop' or param_name == '':
                    Xs_tmp[int(split_line[1])] = int(split_line[2])
                    i += 1

                    l = lines[i]
                    split_line = l.split(',')
                    param_name = split_line[0]
                param_dict['crop'] = Xs_tmp
            elif (param_name == 'time_to_remove' or param_name == 'voxel_downsampled' or
                  param_name == 'original_voxel_ratio' or param_name == 'TP_to_keep'):
                name = param_name
                out = []
                while (name == param_name or param_name == '') and  i < nb_lines:
                    try:
                        out += [int(split_line[1])]
                    except ValueError:
                        out += [float(split_line[1])]
                    i += 1
                    if i < nb_lines:
                        l = lines[i]
                        split_line = l.split(',')
                        param_name = split_line[0]
                param_dict[name] = np.array(out)
        keys_to_have = ['original_voxel_ratio', 'file_name', 'ending_point',
                        'reference_time', 'crop', 'voxel_downsampled',
                        'starting_point', 'folder', 'output_folder']
        assert set(keys_to_have).difference(set(param_dict.keys())) == set(), 'Problem when parsing the param file'

        if 'registered_output' in param_dict:
            registered_folder += [param_dict['registered_output'] + param_dict['registered_file_name']]
        else:
            registered_folder += [[]]
        paths += [param_dict['folder'] + param_dict['file_name']]
        Xs += [param_dict['crop']]
        reference_times += [param_dict['reference_time']]
        time_list_tmp = range(param_dict['starting_point'], param_dict['ending_point'])
        for t in param_dict.get('time_to_remove', [-1]):
            if t in time_list_tmp:
                time_list_tmp.remove(t)
        time_list += [time_list_tmp]
        VS += [param_dict['voxel_downsampled']]
        OR += [param_dict['original_voxel_ratio']]
        suffixes += [param_dict['output_folder']]
        times += [len(time_list_tmp)]
        registration_methods += [param_dict.get('registration_method', 'blockmatching')]
        comp_trsfs += [param_dict.get('compute_trsf', 'yes') == 'yes']
        app_trsfs += [param_dict.get('apply_trsf', 'no') == 'yes']
        TP_to_keep += [param_dict.get('TP_to_keep', [-1])]
        paths_to_bin += [param_dict.get('path_to_bin', '')]

    return (paths, reference_times, Xs, times, VS, OR, time_list, suffixes,
            registration_methods, registered_folder, comp_trsfs, app_trsfs, TP_to_keep, paths_to_bin)

if __name__ == '__main__':
    # reads the parameters for the run
    (paths, reference_times, Xs, times, VS, OR, time_list,
        folder_names, registration_methods, registered_folders,
        comp_trsfs, app_trsfs, TP_to_keeps, paths_to_bin) = read_param_file()

    # ask for the number of cpu
    nb_cpu = int(raw_input('Please enter #cpus desired (-1 implies all, max cpus = %d; %dGB of RAM): '%(cpu_count(), psutil.virtual_memory()[1]/10**9)))
    if nb_cpu < 1:
        nb_cpu = cpu_count()
    
    for (p, r, X, nb_times, vs, ors,
         tp_list, folder_name, registration_method,
         registered_folder, comp_trsf,
         app_trsf, TP_to_keep, path_to_bin) in zip(paths, reference_times,
                                      Xs, times, VS, OR, time_list,
                                      folder_names, registration_methods,
                                      registered_folders, comp_trsfs, app_trsfs, TP_to_keeps, paths_to_bin):
        if not os.path.exists('.TMP'):
            os.makedirs('.TMP')
        first_TP = tp_list[0]
        trsf_p = folder_name + '_trsfs/'
        if comp_trsf:
            print 'Starting registration computation'
            if registration_method == 'blockmatching':
                run_produce_trsf(p, r, X, nb_times, trsf_p, tp_list, ors, vs = vs, first_TP = first_TP, nb_cpu = nb_cpu)
            elif registration_method == 'phase_correlation':
                run_produce_trsf_simple(p, r, X, nb_times, trsf_p, tp_list, ors, vs = vs, first_TP = first_TP, nb_cpu = nb_cpu)
            print 'Transformation composition'
            compose_trsf(first_TP, r, trsf_p, tp_list)
            compose_trsf(tp_list[-1], r, trsf_p, tp_list)
            print_trsf(trsf_p, folder_name + '_out/', nb_times, r, tp_list, vs, ors, first_TP)

        if registered_folder != [] and app_trsf:
            applyTrsf_from_run(trsf_p, nb_times, r, tp_list, vs, ors, first_TP, registered_folder, nb_cpu, TP_to_keep)
            assemble_projection(registered_folder, tp_list, folder_name, proj = 'xy')
            assemble_projection(registered_folder, tp_list, folder_name, proj = 'yz')
            assemble_projection(registered_folder, tp_list, folder_name, proj = 'xz')
        else:
            print 'Registration and Projection'
            p_im = '.TMP/' + p.split('/')[-1].replace('.klb', '.inr')
            p_trsf = trsf_p + '/t%03d-%03d.txt'
            mapping = [(p_im, p_trsf, t, r, vs) for t in tp_list]

            # build all the MIP for every time point
            pool = Pool()
            out = pool.map(run_shift_proj, mapping)
            pool.close()
            pool.terminate()

            XY_shape = np.array(list(np.max([(ii[0].shape) for ii in out], axis = 0)) + [nb_times+1])
            YZ_shape = np.array(list(np.max([(ii[1].shape) for ii in out], axis = 0)) + [nb_times+1])
            XZ_shape = np.array(list(np.max([(ii[2].shape) for ii in out], axis = 0)) + [nb_times+1])
            XY = np.zeros(XY_shape, dtype = np.uint16)
            YZ = np.zeros(YZ_shape, dtype = np.uint16)
            XZ = np.zeros(XZ_shape, dtype = np.uint16)
            for xy, yz, xz, t in out:
                b1, b2 = (XY_shape[:-1] - xy.shape)/2
                XY[b1:b1+xy.shape[0],b2:b2+xy.shape[1],t] = xy
                b1, b2 = (YZ_shape[:-1] - yz.shape)/2
                YZ[b1:b1+yz.shape[0],b2:b2+yz.shape[1],t] = yz
                b1, b2 = (XZ_shape[:-1] - xz.shape)/2
                XZ[b1:b1+xz.shape[0],b2:b2+xz.shape[1],t] = xz

            imsave(folder_name + '_out/XY_proj.tif', SpatialImage(XY))
            imsave(folder_name + '_out/YZ_proj.tif', SpatialImage(YZ))
            imsave(folder_name + '_out/XZ_proj.tif', SpatialImage(XZ))

    shutil.rmtree('.TMP')