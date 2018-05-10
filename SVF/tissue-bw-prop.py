# This file is subject to the terms and conditions defined in
# file 'LICENSE', which is part of this source code package.
# Author: Leo Guignard (guignardl...@AT@...janelia.hhmi.org)

from IO import imread, imsave, SpatialImage
from scipy import ndimage as nd
import numpy as np
import os
from multiprocessing import Pool
from TGMMlibraries import lineageTree
from scipy import interpolate
import sys

def get_spherical_coordinates(x, y, z):
    ''' Computes spherical coordinates for an x, y, z Cartesian position
    '''
    r = np.linalg.norm([x, y, z])
    theta = np.arctan2(y, x)
    phi = np.arccos(z/r)
    alpha = (np.pi/2 + np.arctan2(x, z)) % (2*np.pi)
    return r, theta, phi, alpha

def write_header_am_2(f, nb_points, length):
    ''' Header for Amira .am files
    '''
    f.write('# AmiraMesh 3D ASCII 2.0\n')
    f.write('define VERTEX %d\n'%(nb_points*2))
    f.write('define EDGE %d\n'%nb_points)
    f.write('define POINT %d\n'%((length)*nb_points))
    f.write('Parameters {\n')
    f.write('\tContentType "HxSpatialGraph"\n')
    f.write('}\n')

    f.write('VERTEX { float[3] VertexCoordinates } @1\n')
    f.write('EDGE { int[2] EdgeConnectivity } @2\n')
    f.write('EDGE { int NumEdgePoints } @3\n')
    f.write('POINT { float[3] EdgePointCoordinates } @4\n')
    f.write('VERTEX { float Vcolor } @5\n')
    f.write('VERTEX { int Vbool } @6\n')
    f.write('EDGE { float Ecolor } @7\n')
    f.write('VERTEX { int Vbool2 } @8\n')

def write_to_am_2(path_format, LT_to_print, t_b = None, t_e = None, length = 5, manual_labels = None, 
                  default_label = 5, new_pos = None):
    ''' Writes a lineageTree into an Amira readable data (.am format).
        Args:
            path_format: string, path to the output. It should contain 1 %03d where the time step will be entered
            LT_to_print: lineageTree, lineageTree to write
            t_b: int, first time point to write (if None, min(LT.to_take_time) is taken)
            t_e: int, last time point to write (if None, max(LT.to_take_time) is taken)
                note: if there is no 'to_take_time' attribute, LT_to_print.time_nodes is considered instead
                    (historical)
            length: int, length of the track to print (how many time before).
            manual_labels: {id: label, }, dictionary that maps cell ids to 
            default_label: int, default value for the manual label
            new_pos: {id: [x, y, z]}, dictionary that maps a 3D position to a cell ID.
                if new_pos == None (default) then LT_to_print.pos is considered.
    '''
    if not hasattr(LT_to_print, 'to_take_time'):
        LT_to_print.to_take_time = LT_to_print.time_nodes
    if t_b is None:
        t_b = min(LT_to_print.to_take_time.keys())
    if t_e is None:
        t_e = max(LT_to_print.to_take_time.keys())
    if new_pos is None:
        new_pos = LT_to_print.pos

    if manual_labels is None:
        manual_labels = {}
    for t in range(t_b, t_e + 1):
        f = open(path_format%t, 'w')
        nb_points = len(LT_to_print.to_take_time[t])
        write_header_am_2(f, nb_points, length)
        points_v = {}
        for C in LT_to_print.to_take_time[t]:
            C_tmp = C
            positions = []
            for i in xrange(length):
                C_tmp = LT_to_print.predecessor.get(C_tmp, [C_tmp])[0]
                positions.append(new_pos[C_tmp])
            points_v[C] = positions

        f.write('@1\n')
        for i, C in enumerate(LT_to_print.to_take_time[t]):
            f.write('%f %f %f\n'%tuple(points_v[C][0]))
            f.write('%f %f %f\n'%tuple(points_v[C][-1]))

        f.write('@2\n')
        for i, C in enumerate(LT_to_print.to_take_time[t]):
            f.write('%d %d\n'%(2*i, 2*i+1))

        f.write('@3\n')
        for i, C in enumerate(LT_to_print.to_take_time[t]):
            f.write('%d\n'%(length))

        f.write('@4\n')
        tmp_velocity = {}
        for i, C in enumerate(LT_to_print.to_take_time[t]):
            for p in points_v[C]:
                f.write('%f %f %f\n'%tuple(p))

        f.write('@5\n')
        for i, C in enumerate(LT_to_print.to_take_time[t]):
            f.write('%f\n'%(manual_labels.get(C, default_label)))
            f.write('%f\n'%(0))

        f.write('@6\n')
        for i, C in enumerate(LT_to_print.to_take_time[t]):
            f.write('%d\n'%(int(manual_labels.get(C, default_label) != default_label)))
            f.write('%d\n'%(0))
        
        f.write('@7\n')
        for i, C in enumerate(LT_to_print.to_take_time[t]):
            f.write('%f\n'%(np.linalg.norm(points_v[C][0] - points_v[C][-1])))

        f.write('@8\n')
        for i, C in enumerate(LT_to_print.to_take_time[t]):
            f.write('%d\n'%(1))
            f.write('%d\n'%(0))
        f.close()

def read_param_file():
    ''' Asks for, reads and formats the parameter file
    '''
    p_param = raw_input('Please enter the path to the parameter file/folder:\n')
    p_param = p_param.replace('"', '')
    p_param = p_param.replace("'", '')
    p_param = p_param.replace(" ", '')
    if p_param[-4:] == '.csv':
        f_names = [p_param]
    else:
        f_names = [os.path.join(p_param, f) for f in os.listdir(p_param) if '.csv' in f and not '~' in f]
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
            if param_name in ['labels', 'downsampling']:
                name = param_name
                out = []
                while (name == param_name or param_name == '') and  i < nb_lines:
                    if split_line[1].isdigit():
                        out += [int(split_line[1])]
                    else:
                        out += [float(split_line[1])]
                    i += 1
                    if i < nb_lines:
                        l = lines[i]
                        split_line = l.split(',')
                        param_name = split_line[0]
                param_dict[name] = np.array(out)
            elif param_name in ['label_names']:
                name = param_name
                out = []
                while (name == param_name or param_name == '') and  i < nb_lines:
                    out += [split_line[1].replace('\n', '').replace('\r', '')]
                    i += 1
                    if i < nb_lines:
                        l = lines[i]
                        split_line = l.split(',')
                        param_name = split_line[0]
                param_dict[name] = np.array(out)
            else:
                param_dict[param_name] = split_line[1].strip()
                i += 1
            if param_name == 'time':
                param_dict[param_name] = int(split_line[1])
        path_LT = param_dict.get('path_to_LT', '.')
        path_VF = param_dict.get('path_to_VF', '.')
        path_mask = param_dict.get('path_to_mask', '.')
        t = param_dict.get('time', 0)
        path_out_am = param_dict.get('path_to_am', '.')
        labels = param_dict.get('labels', [])
        DS = param_dict.get('downsampling', [])
        ani = np.float(param_dict.get('anisotropy', 1.))
        path_DB = param_dict.get('path_DB', '.')
        path_div = param_dict.get('path_div', None)
        path_bary = param_dict.get('path_bary', None)
        label_names = param_dict.get('label_names', None)
        invert = param_dict.get('inverted', '1') != '0'

    return (path_LT, path_VF, path_mask, t, path_out_am,
            labels, DS, path_DB, path_div, path_bary,
            label_names, ani, invert)

def get_division_mapping(path_div, VF):
    ''' Computes the mapping between found divisions and SVF objects
        Args:
            path_div: sting, name of the division file
            VF: lineageTree
    '''
    ass_div = {}
    if path_div is not None:
        f = open(path_div)
        lines = f.readlines()
        f.close()
        divisions_per_time = {}
        for l in lines[1:]:
            x, y, z, t = np.array(l.split(',')[:-1]).astype(float)
            if t in VF.time_nodes:
                divisions_per_time.setdefault(int(t), []).append(np.array([x, y, z]) * [1, 1, 5])

        div_in_VF = {}
        dist_to_div = {}
        for t, d in divisions_per_time.iteritems():
            idx3d, data = VF.get_idx3d(t)
            dist, idxs = idx3d.query(d)
            div_C = np.array(data)[idxs]
            dist_to_div.update(dict(zip(div_C, dist)))
            ass_div.update(dict(zip(div_C, d)))
    return ass_div

def write_DB(path_DB, path_div, VF, tracking_value, tb, te):
    ''' Write the csv database in Database.csv
        Args:
            path_DB: string, path to the output database
            path_div: string, path to the potential division file
            VF: lineageTree
            tracking_value: {int: int, }, dictionary that maps an object id to a label
            tb: int, first time point to write
            te: int, last time point to write
    '''
    ass_div = get_division_mapping(path_div, VF)
    f2 = open(path_DB + 'Database-Tissues.csv', 'w')
    f2.write('id, mother_id, x, y, z, r, theta, phi, t, label, D-x, D-y, D-z, D-r, D-theta, D-phi\n')
    for t in range(tb, te+1):
        for c in VF.time_nodes[t]:
            S_p = (-1, -1, -1)
            if VF.predecessor.get(c, []) != []:
                M_id = VF.predecessor[c][0]
            else:
                M_id = -1
            P = tuple(VF.pos[c])
            if path_bary is not None:
                S_p = tuple(get_spherical_coordinates(*(barycenters[t] - VF.pos[c]))[:-1])
            L = tracking_value.get(c, -1)
            D_P = tuple(ass_div.get(c, [-1, -1, -1]))
            if path_bary is not None:
                D_S_p = (-1, -1, -1) if not c in ass_div else tuple(get_spherical_coordinates(*(barycenters[t] - ass_div[c]))[:-1])
            else:
                D_S_p = (-1, -1, -1)
            f2.write(('%d, %d, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %d, %d,' + 
                      '%.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n')%((c, M_id) + P + S_p + (t, L) + D_P + D_S_p))
    f2.close()

def get_barycenter(fname, tb, te):
    ''' Reads and coes a linear piecewise interpolation/extrapolation barycenters
        Args:
            fname: string, name of the barycenter file (each line as 'x, y, z, t')
            tb: first time point to interpolate
            te: last time point to interpolate
        Returns:
            barycenters_interp: {int:[float, float, float], }, dictionary mapping
                        a time point to the interpolated barycenter at that time
            barycenters: {int:[float, float, float], }, dictionary mapping
                        a time point to the barycenter for each time in fname
    '''
    f = open(fname)
    lines = f.readlines()
    f.close()
    barycenters = {}
    for l in lines:
        split_l = l.split(',')
        try:
            barycenters[int(split_l[-1])] = tuple(float(v) for v in split_l[:-1])
        except Exception as e:
            pass
    times = sorted(barycenters)
    Xb, Yb, Zb = np.array([barycenters[t] for t in times]).T
    Xb_f = interpolate.InterpolatedUnivariateSpline(times, Xb, k=1)
    Yb_f = interpolate.InterpolatedUnivariateSpline(times, Yb, k=1)
    Zb_f = interpolate.InterpolatedUnivariateSpline(times, Zb, k=1)
    Ti = np.arange(tb - 1, te + 2)
    barycenters_interp = dict(zip(Ti, zip(Xb_f(Ti), Yb_f(Ti), Zb_f(Ti))))

    return barycenters_interp, barycenters


if __name__ == '__main__':
    (path_LT, path_VF, path_mask, t, path_out_am,
     labels, DS, path_DB, path_div, path_bary,
     label_names, ani, invert) = read_param_file()
    if not os.path.exists(path_out_am):
        os.makedirs(path_out_am)
    if not os.path.exists('.mask_images/'):
        os.makedirs('.mask_images/')
    VF = lineageTree(path_VF)
    tb = VF.t_b
    te = VF.t_e

    if path_bary is not None:
        try:
            barycenters, b_dict = get_barycenter(path_bary, tb, te)
        except Exception as e:
            print "Wrong file path to barycenter, please specify the path to the .csv file."
            print "The process will continue as if no barycenter were provided,"
            print "disabling the computation of the spherical coordinates"
            print "error raised: ", e
            path_bary = None

    im = imread(path_mask)
    for l in labels:
        masked_im = im == l
        tmp = nd.binary_opening(masked_im, iterations = 3)
        tmp = nd.binary_closing(tmp, iterations = 4)
        imsave('.mask_images/%03d.tif'%l, SpatialImage(tmp).astype(np.uint8))

    mask_dir = '.mask_images/'
    if label_names is not None:
        masks = sorted([('.mask_images/%03d.tif'%l, label_names[i]) for i, l in enumerate(labels)], cmp=lambda x1, x2:cmp(x1[1], x2[1]))
        masks = [m[0] for m in masks]
    else:
        masks = ['.mask_images/%03d.tif'%l for l in labels] 

    init_cells = {m: set() for m in range(len(masks))}
    x_max, y_max, z_max = 0, 0, 0

    for i, path_mask in enumerate(masks):
        if invert:
            mask = imread(path_mask).transpose(1, 0, 2)
            mask = mask[:,::-1,:]
        else:
            mask = imread(path_mask)
        max_vals = np.array(mask.shape) - 1
        for c in VF.time_nodes[t]:
            pos_rounded = np.floor(VF.pos[c]/(np.array(DS)*[1.,1.,ani])).astype(np.int)
            pos_rounded = tuple(np.min([max_vals, pos_rounded], axis = 0))
            if mask[pos_rounded]:
                init_cells[i].add(c)

    tracking_value = {}
    for t, cs in init_cells.iteritems():
        for c in cs:
            to_treat = [c]
            tracking_value.setdefault(c, set()).add(t)
            while to_treat != []:
                c_tmp = to_treat.pop()
                next_cells = VF.successor.get(c_tmp, [])
                to_treat += next_cells
                for n in next_cells:
                    tracking_value.setdefault(n, set()).add(t)
            to_treat = [c]
            tracking_value.setdefault(c, set()).add(t)
            while to_treat != []:
                c_tmp = to_treat.pop()
                next_cells = VF.predecessor.get(c_tmp, [])
                to_treat += next_cells
                for n in next_cells:
                    tracking_value.setdefault(n, set()).add(t)

    tracking_value = {k:np.sum(list(v)) for k, v in tracking_value.iteritems() if len(v) == 1}

    write_to_am_2(path_out_am + '/seg_t%04d.am', VF, t_b= tb, t_e= te,
                manual_labels = tracking_value, default_label = np.max(tracking_value.values())+1,
                length = 7)

    for im_p in masks:
        os.remove(im_p)

    write_DB(path_DB, path_div, VF, tracking_value, tb, te)
