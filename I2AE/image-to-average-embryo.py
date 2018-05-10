#!python
# This file is subject to the terms and conditions defined in
# file 'LICENCE', which is part of this source code package.
# Author: Leo Guignard (guignardl...@AT@...janelia.hhmi.org)

import scipy as sp
from scipy import interpolate
from itertools import product
import re
import sys
import os
from TGMMlibraries import lineageTree
from time import time
import struct
from multiprocessing import Pool
from itertools import combinations
# import xml.etree.ElementTree as ET
import numpy as np
from scipy.spatial import cKDTree as kdtree
from scipy.spatial import Delaunay
from scipy import spatial
import itertools as it
from scipy import spatial
from scipy.optimize import linear_sum_assignment
from copy import deepcopy
from scipy import ndimage as nd
import cPickle as pkl
from IO import imread, imsave, SpatialImage

def rigid_transform_3D(A, B):
    ''' Compute the 4x4 matrix reprenting the 3D rigid rotation
        minimizing the least squares between two paired sets of points *A* and *B*
        Args:
            A: Nx3 array, the first set of 3D points
            B: Nx3 array, the second set of 3D points
        Returns:
            out: 4x4 array, 3D rigid rotation
    '''
    assert len(A) == len(B)
    
    if not type(A) == np.matrix:
        A = np.matrix(A)

    if not type(B) == np.matrix:
        B = np.matrix(B)

    N = A.shape[0]; # total points

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)
    
    # centre the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = np.transpose(AA) * BB

    U, S, Vt = np.linalg.svd(H)

    R = Vt.T * U.T

    # special reflection case
    if np.linalg.det(R) < 0:
        Vt[2,:] *= -1
        R = Vt.T * U.T

    t = -R*centroid_A.T + centroid_B.T

    out = np.identity(4)
    out[:3, :3] = R
    out[:-1, 3:] = t
    return out


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
                  default_label = 5, new_pos = None, to_take_time = None):
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
            to_take_time: {t, [int, ]}
    '''
    if to_take_time is None:
        to_take_time = LT_to_print.time_nodes
    if t_b is None:
        t_b = min(to_take_time.keys())
    if t_e is None:
        t_e = max(to_take_time.keys())
    if new_pos is None:
        new_pos = LT_to_print.pos

    if manual_labels is None:
        manual_labels = {}
    for t in range(t_b, t_e + 1):
        f = open(path_format%t, 'w')
        nb_points = len(to_take_time[t])
        write_header_am_2(f, nb_points, length)
        points_v = {}
        for C in to_take_time[t]:
            C_tmp = C
            positions = []
            for i in xrange(length):
                c_before = C_tmp
                C_tmp = LT_to_print.predecessor.get(C_tmp, [C_tmp])[0]
                while not C_tmp in new_pos and C_tmp in LT_to_print.predecessor:
                    C_tmp = LT_to_print.predecessor.get(C_tmp, [C_tmp])[0]
                if not C_tmp in new_pos:
                    C_tmp = c_before
                positions.append(np.array(new_pos[C_tmp]))
            points_v[C] = positions

        f.write('@1\n')
        for i, C in enumerate(to_take_time[t]):
            f.write('%f %f %f\n'%tuple(points_v[C][0]))
            f.write('%f %f %f\n'%tuple(points_v[C][-1]))

        f.write('@2\n')
        for i, C in enumerate(to_take_time[t]):
            f.write('%d %d\n'%(2*i, 2*i+1))

        f.write('@3\n')
        for i, C in enumerate(to_take_time[t]):
            f.write('%d\n'%(length))

        f.write('@4\n')
        tmp_velocity = {}
        for i, C in enumerate(to_take_time[t]):
            for p in points_v[C]:
                f.write('%f %f %f\n'%tuple(p))

        f.write('@5\n')
        for i, C in enumerate(to_take_time[t]):
            f.write('%f\n'%(manual_labels.get(C, default_label)))
            f.write('%f\n'%(0))

        f.write('@6\n')
        for i, C in enumerate(to_take_time[t]):
            f.write('%d\n'%(int(manual_labels.get(C, default_label) != default_label)))
            f.write('%d\n'%(0))
        
        f.write('@7\n')
        for i, C in enumerate(to_take_time[t]):
            f.write('%f\n'%(np.linalg.norm(points_v[C][0] - points_v[C][-1])))

        f.write('@8\n')
        for i, C in enumerate(to_take_time[t]):
            f.write('%d\n'%(1))
            f.write('%d\n'%(0))

        f.close()

def write_to_am_rem(path_format, LT_to_print, t_b = None, t_e = None, length = 5, manual_labels = None, 
                  default_label = 5, new_pos = None, to_take_time = None, to_remove = None, predecessor = None):
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
            to_take_time: {t, [int, ]}
            to_remove: [int, ], list of id not to take into account
            predecessor: {int: [int, ]}, dictionary that maps a cell id to the list of its
                predecesors
    '''
    if to_take_time is None:
        to_take_time = LT_to_print.time_nodes
    if t_b is None:
        t_b = min(to_take_time.keys())
    if t_e is None:
        t_e = max(to_take_time.keys())
    if new_pos is None:
        new_pos = LT_to_print.pos
    if to_remove is None:
        to_remove = set()
    if predecessor is None:
        predecessor = LT_to_print.predecessor

    if manual_labels is None:
        manual_labels = {}
    for t in range(t_b, t_e + 1):
        f = open(path_format%t, 'w')
        Points = [c for c in to_take_time[t] if not c in to_remove]
        nb_points = len(Points)
        write_header_am_2(f, nb_points, length)
        points_v = {}
        for C in Points:
            C_tmp = C
            positions = []
            for i in xrange(length):
                c_before = C_tmp
                C_tmp = predecessor.get(C_tmp, [C_tmp])[0]
                while not C_tmp in new_pos and C_tmp in predecessor:
                    C_tmp = predecessor.get(C_tmp, [C_tmp])[0]
                if not C_tmp in new_pos:
                    C_tmp = c_before
                positions.append(np.array(new_pos[C_tmp]))
            points_v[C] = positions

        f.write('@1\n')
        for i, C in enumerate(Points):
            f.write('%f %f %f\n'%tuple(points_v[C][0]))
            f.write('%f %f %f\n'%tuple(points_v[C][-1]))

        f.write('@2\n')
        for i, C in enumerate(Points):
            f.write('%d %d\n'%(2*i, 2*i+1))

        f.write('@3\n')
        for i, C in enumerate(Points):
            f.write('%d\n'%(length))

        f.write('@4\n')
        tmp_velocity = {}
        for i, C in enumerate(Points):
            for p in points_v[C]:
                f.write('%f %f %f\n'%tuple(p))

        f.write('@5\n')
        for i, C in enumerate(Points):
            f.write('%f\n'%(manual_labels.get(C, default_label)))
            f.write('%f\n'%(0))

        f.write('@6\n')
        for i, C in enumerate(Points):
            f.write('%d\n'%(int(manual_labels.get(C, default_label) != default_label)))
            f.write('%d\n'%(0))
        
        f.write('@7\n')
        for i, C in enumerate(Points):
            f.write('%f\n'%(np.linalg.norm(points_v[C][0] - points_v[C][1])))
            # if 15 < np.linalg.norm(points_v[C][0] - points_v[C][1]):
            #     print C

        f.write('@8\n')
        for i, C in enumerate(Points):
            f.write('%d\n'%(1))
            f.write('%d\n'%(0))

        f.close()

def applyTrsf_full(mats, p, t):
    ''' Provided a dictionary of matrices, a position and a time,
        applies the appropriate transformation to the position *p*
        Args:
            mats: {int: 4x4 array, }, dictionary where the key is time
                and the value is the associated transformation for that time
            p: [float, float, float], x, y, z coordinate in the Cartesian referential
            t: int, time to which p belongs to
        Returns:
            new_p: [float, float, float], the new position of the original point *p*
    '''
    if not t in mats:
        t = mats.keys()[np.argmin(np.abs(t - np.array(mats.keys())))]
    new_p = applyMatrixTrsf(p, mats[t])
    return new_p

def applyTrsf_time_laps(LT, timeTrsf, mats):
    ''' Apply the appropriate time and rigid Transformations
        to the points in a lineageTree provided the time and rigid transformations.
        If necessary, time points are added/removed.
        The nodes are stored in LT.new_time_nodes for there right times,
        The transformed positions are stored in LT.new_pos
        Args:
            LT: lineageTree
            timeTrsf: time transformation function
            mats: {int: 4x4 array, }, dictionary where the key is time
                and the value is the associated transformation for that time
    '''
    LT.new_time_nodes = {}
    LT.new_pos = {}
    t_total = []
    t_to_fill = []
    time_evol = []
    if not hasattr(LT, 'to_take_time'):
        LT.to_take_time = LT.time_nodes

    time_mapping = np.round(timeTrsf(sorted(LT.time_nodes.keys()))).astype(np.int)
    times_to_fuse = time_mapping[1:][(time_mapping[1:] == time_mapping[:-1])]
    fill = zip(time_mapping[:-1][1 < (time_mapping[1:] - time_mapping[:-1])],
               time_mapping[1:][1 < (time_mapping[1:] - time_mapping[:-1])])

    cells_removed = set()
    treated = set()
    for t in sorted(LT.time_nodes.keys()):
        if not t in treated:
            treated.add(t)
            new_t = int(np.round(timeTrsf(t)))
            LT.new_time_nodes[new_t] = []
            if not new_t in times_to_fuse:
                for c in LT.to_take_time[t]:
                    LT.new_pos[c] = np.array(applyTrsf_full(mats, LT.pos[c], t))
                    LT.new_time_nodes[new_t].append(c)
            else:
                ti = t
                t_to_aggregate = []
                while int(np.round(timeTrsf(ti))) == new_t:
                    t_to_aggregate += [ti]
                    treated.add(ti)
                    ti += 1
                for ci in LT.to_take_time[t_to_aggregate[0]]:
                    track = [ci]
                    while LT.successor.get(track[-1], []) != [] and LT.time[track[-1]] < t_to_aggregate[-1]:
                        track += LT.successor[track[-1]]
                    glob_pred = -1
                    glob_pred = -1
                    if (LT.time.get(LT.successor.get(track[-1], [-1])[0], -1) == t_to_aggregate[-1] + 1 and
                        LT.predecessor.get(track[0], [-1])[0] != -1):
                        all_pos = [applyTrsf_full(mats, LT.pos[cii], LT.time[cii]) for cii in track]
                        avg_p = np.mean(all_pos, axis = 0)
                        glob_pred = LT.predecessor[track[0]][0]
                        glob_succ = LT.successor[track[-1]][0]
                    for cii in track:
                        cells_removed.add(cii)
                        cii_succ = LT.successor.pop(cii, [])
                        for succ in cii_succ:
                            if cii in LT.predecessor.get(succ, []):
                                LT.predecessor[succ].remove(cii)
                            if LT.predecessor.get(succ, []) == []:
                                LT.predecessor.pop(succ)
                        cii_pred = LT.predecessor.pop(cii, [])
                        for pred in cii_pred:
                            if cii in LT.successor.get(pred, []):
                                LT.successor[pred].remove(cii)
                            if LT.successor.get(pred, []) == []:
                                LT.successor.pop(pred)
                    if glob_pred != -1 and glob_succ != -1:
                        C_id = LT.get_next_id()
                        LT.successor[pred] = [C_id]
                        LT.predecessor[C_id] = [glob_pred]
                        LT.successor[C_id] = [glob_succ]
                        LT.predecessor[succ] = [C_id]
                        LT.edges.append((C_id, glob_succ))
                        LT.edges.append((glob_pred, C_id))
                        LT.nodes.append(C_id)
                        LT.new_time_nodes[new_t] += [C_id]
                        LT.new_pos[C_id] = avg_p
    for tb, te in fill:
        for t_inter in range(tb + 1, te):
            LT.new_time_nodes[t_inter] = []
        for ci in LT.new_time_nodes[tb]:
            next_C = LT.successor.get(ci, [-1])[0]
            if next_C != -1:
                p_init = LT.new_pos[ci]
                p_final = LT.new_pos[LT.successor[ci][0]]
                vals = np.zeros((2, 4))
                vals[0, :-1] = p_init
                vals[0, -1] = tb
                vals[1, :-1] = p_final
                vals[1, -1] = te
                pos_interp = interp_3d_coord(vals)
                c_pred = ci
                for t_inter in range(tb + 1, te):
                    C_id = LT.get_next_id()
                    LT.nodes += [C_id]
                    LT.new_time_nodes[t_inter] += [C_id]
                    LT.successor[c_pred] = [C_id]
                    LT.predecessor[C_id] = [c_pred]
                    LT.edges.append((c_pred, C_id))
                    c_pred = C_id
                    LT.new_pos[C_id] = pos_interp(t_inter)
                LT.successor[c_pred] = [next_C]
                LT.predecessor[next_C] = [c_pred]
                LT.edges.append((c_pred, C_id))
                # LT.edges.remove((ci, next_C))
    LT.nodes = list(set(LT.nodes).difference(cells_removed))


def get_points_from_ray(points, A, B, r = 5):
    ''' Build the projections of points on a line AB that are at a distance
        maximum of *r*. It return three lists.
        Args:
            points: [[float, float, float], ], list of 3D positions
            A: [float, float, float], 3D position, first point of the line
            B: [float, float, float], 3D position, second point of the line
            r: float, distance maximum between the line and the points
        Returns:
            [[float, float, float], ], list for position of the selected points after projection
            [float, ], list of distances of the selected points to *A*
            [[float, float, float], ], list for position of the selected points before projection
    '''
    AB = B - A
    Ap = points - A
    n_Aproj = (np.dot(Ap, AB)/np.linalg.norm(AB)).reshape(len(Ap), 1)
    proj = (A + n_Aproj * AB/np.linalg.norm(AB))
    distances = np.linalg.norm(points - proj, axis = 1)
    distr_on_AB = n_Aproj.reshape(len(Ap))
    return (proj[(distances<r) & (distr_on_AB >= 0)],
            distr_on_AB[(distances<r) & (distr_on_AB >= 0)],
            points[(distances<r) & (distr_on_AB >= 0)])

N_func = lambda x, mu, sigma: np.exp(-(x - mu)**2/(2*sigma**2)) / (np.sqrt(2*np.pi*sigma**2))
def get_low_and_high(points, A, B, r = 150, percentile = 80):
    ''' Retrieves the high and low points along a ray casted from a point A to a point B
        Args:
            points: [[float, float, float], ], list of 3D positions
            A: [float, float, float], 3D position, first point of the line
            B: [float, float, float], 3D position, second point of the line
            r: float, distance maximum between the line and the points
            percentile: float, percentile of the distribution kept, has to be between 0 and 100.
    '''
    kept_points, distribution_on_AB, kept_full_points = get_points_from_ray(points, A, B, 150)
    m, s = np.mean(distribution_on_AB), np.std(distribution_on_AB)
    mask = ((m - 5 * s) < distribution_on_AB) & (distribution_on_AB < (m + 5 * s))
    distribution_on_AB = distribution_on_AB[mask]
    kept_points = kept_points[mask]
    if len(distribution_on_AB)>10:
        th_v_m = np.percentile(distribution_on_AB, percentile)
        th_v_p = np.percentile(distribution_on_AB, 100-percentile)

        low_p = kept_points[np.argmin(np.abs(distribution_on_AB - th_v_m))]
        high_p = kept_points[np.argmin(np.abs(distribution_on_AB - th_v_p))]
        return low_p, high_p, distribution_on_AB, kept_points, kept_full_points
    else:
        return None, None,None, None, None

percentile = 85
def build_points(t, A):
    ''' Builds the inner and outer shell of the two lineage trees LT_1 and LT_2
        for time point *t*t from rays casted from the center of mass *A*.
        The points are written in the files:
        pts01_high_DS_%03d.txt, pts01_low_DS_%03d.txt, pts02_high_DS_%03d.txt, pts02_low_DS_%03d.txt
        Args:
            t: int, time point to treat
            A: [float, float, float], 3D coordinates of the center of mass
    '''
    tic = time()
    A = np.array(A)
    points_2 = np.array(new_pos_avg_TGMM_rot.values())
    points_1 = np.array([ref_TGMM.pos[c] for c in ref_TGMM.time_nodes[t]])
    low_points_2 = []
    high_points_2 = []
    low_points_1 = []
    high_points_1 = []
    tmp_pts = []
    tmp_dist = []
    nb_p_for_phi = lambda phi: np.round(np.interp(phi, [0, np.pi/2], [1, 50])).astype(np.int)
    for phi in np.linspace(0, np.pi/2, 50):
        for theta in np.linspace(0, np.pi/2, nb_p_for_phi(phi)):
            Z = 500*np.cos(theta)*np.sin(phi)
            Y = 500*np.sin(theta)*np.sin(phi)
            X = 500*np.cos(phi)
            sphere_points = [[ X,  Y,  Z], [ X,  Y, -Z], [ X, -Y,  Z], [ X, -Y, -Z]]
            for B in sphere_points:
                tmp_pts += [A - B]
                low_p_2, high_p_2, distribution_on_AB_2, kept_points_2, kept_full_points_2 = get_low_and_high(points_2, A, A - B, 150, percentile)
                low_p_1, high_p_1, distribution_on_AB_1, kept_points_1, kept_full_points_1 = get_low_and_high(points_1, A, A - B, 150, percentile)
                if (not low_p_1 is None) and (not low_p_2 is None):
                    tmp_dist += [[distribution_on_AB_2, kept_points_2, kept_full_points_2]]
                    low_points_2 += [low_p_2]
                    high_points_2 += [high_p_2]
                    low_points_1 += [low_p_1]
                    high_points_1 += [high_p_1]
                low_p_2, high_p_2, distribution_on_AB_2, kept_points_2, kept_full_points_2 = get_low_and_high(points_2, A, A + B, 150, percentile)
                low_p_1, high_p_1, distribution_on_AB_1, kept_points_1, kept_full_points_1 = get_low_and_high(points_1, A, A + B, 150, percentile)
                if (not low_p_1 is None) and (not low_p_2 is None):
                    low_points_2.append(low_p_2)
                    high_points_2.append(high_p_2)
                    low_points_1.append(low_p_1)
                    high_points_1.append(high_p_1)

    np.savetxt(match_points_folder + 'pts01_high_DS_%03d.txt'%t, np.array(high_points_1)/DS_value)
    np.savetxt(match_points_folder + 'pts01_low_DS_%03d.txt'%t, np.array(low_points_1)/DS_value)
    np.savetxt(match_points_folder + 'pts02_high_DS_%03d.txt'%t, np.array(high_points_2)/DS_value)
    np.savetxt(match_points_folder + 'pts02_low_DS_%03d.txt'%t, np.array(low_points_2)/DS_value)
    print 't%03d: %.2f'%(t, time() - tic)

def interp_3d_coord(vals, ext = 3):
    ''' Does a piecewise lineare interpolation between 3D points
        Args:
            vals: [[float, float, float, int], ],
                list of 3D values plus there associatedf time
    '''
    if len(vals) > 1:
        vals = vals[np.argsort(vals[:,-1])]
        X_som_front = sp.interpolate.InterpolatedUnivariateSpline(vals[:, -1], vals[:, 0], k=1, ext = ext)
        Y_som_front = sp.interpolate.InterpolatedUnivariateSpline(vals[:, -1], vals[:, 1], k=1, ext = ext)
        Z_som_front = sp.interpolate.InterpolatedUnivariateSpline(vals[:, -1], vals[:, 2], k=1, ext = ext)

        return lambda t: np.array([X_som_front(t), Y_som_front(t), Z_som_front(t)])
    elif ext == 3:
        return lambda t: deepcopy(vals[0][:-1])
    else:
        return lambda t: deepcopy(vals[0][:-1]) if t == vals[0][-1] else np.array([0., 0., 0.])

def get_somite_pos(coord):
    ''' Extract the average position between left and right somites
        Args:
            coord: {string: [float, float, float], }, dictionary that maps the name
                of a landmark to its position
    '''
    s_end_pairing = {}
    for name, p in coord.iteritems():
        if 'LS' in name:
            name_R = name.replace('L', 'R')
            if name_R in coord:
                s_end_pairing.setdefault(int(name[2]), []).append([p, coord[name_R]])
    s_end_pairing = {t:np.mean(v, axis = 1) for t, v in s_end_pairing.iteritems()}
    final_somites = {}
    for S_id in sorted(s_end_pairing):
        final_somites[S_id] = {}
        for S in s_end_pairing[S_id]:
            final_somites[S_id][np.round(S[-1])] = S[:-1]
    return final_somites


def get_somite_front(coord):
    ''' Extract the average position between left and right
        somites at the first time they appear
        Args:
            coord: {string: [float, float, float], }, dictionary that maps the name
                of a landmark to its position
    '''
    s1_paring = []
    for name, p in coord.iteritems():
        if 'LS1.' in name:
            name_R = name.replace('L', 'R')
            if name_R in coord:
                s1_paring.append((p, coord[name_R]))

    somite_front = np.mean(s1_paring, axis = 1)
    return interp_3d_coord(somite_front, ext = 1), {k[-1]: k[:-1] for k in somite_front}

def get_somite_pairing(coord):
    ''' Extract the average position between left and right somites
        Args:
            coord: {string: [float, float, float], }, dictionary that maps the name
                of a landmark to its position
    '''
    corres_N = {'HFPP': 10,
                'AMN': 11,
                'ALM': 12,}
    s_end_pairing = {}
    for name, p in coord.iteritems():
        if 'LS' in name:
            name_R = name.replace('L', 'R')
            if name_R in coord:
                s_end_pairing.setdefault(int(name[2]), []).append([p, coord[name_R]])
        # elif 'HFPP' in name or 'AMN' in name or 'ALM' in name:
        #     s_end_pairing.setdefault(corres_N[name], []).append([p])
    s_end_pairing_mean = {t:np.mean(v, axis = 1) for t, v in s_end_pairing.iteritems()}
    return (np.array([ki[:-1] for k in s_end_pairing_mean.values() for ki in k]),
            sum(s_end_pairing.values(), []))

def get_somite_end(coord):
    ''' Extract the average position between left and right somites
        For the last time the somites have been marked
        Args:
            coord: {string: [float, float, float], }, dictionary that maps the name
                of a landmark to its position
    '''
    s_end_pairing = {}
    for name, p in coord.iteritems():
        if 'LS' in name:
            name_R = name.replace('L', 'R')
            if name_R in coord:
                s_end_pairing.setdefault(int(name[2]), []).append([p, coord[name_R]])
    s_end_pairing = {t:np.mean(v, axis = 1) for t, v in s_end_pairing.iteritems()}
    final = {}
    s_start = {np.int(np.round(np.min(v[:,-1]))):t for t, v in s_end_pairing.iteritems()}
    for s, vals in s_end_pairing.iteritems():
        final[s] = interp_3d_coord(vals)
    choosed_somite = {}
    times = sorted(s_start)
    end_time = times[-1]
    for start in times:
        for t in range(start, end_time + 1):
            choosed_somite.setdefault(t, []).append(s_start[start])

    return final, choosed_somite

def get_somite_interp(coord, somite_to_look):
    s_end_pairing = []
    for name, p in coord.iteritems():
        if 'LS' in name:
            name_R = name.replace('L', 'R')
            if name_R in coord:
                if somite_to_look[np.round((p[-1] + coord[name_R][-1])/2.)] == int(name[2]):
                    s_end_pairing += [np.mean([p, coord[name_R]], axis = 0)]

    return interp_3d_coord(np.array(s_end_pairing), ext = 1), {k[-1]: k[:-1] for k in s_end_pairing}

def get_common_somites_end(coord_1, coord_2_trsf):
    ''' Extract the average position between left and right somites
        Args:
            coord: {string: [float, float, float], }, dictionary that maps the name
                of a landmark to its position
    '''
    somite_end_1, somite_present_1 = get_somite_end(coord_1)
    somite_end_2, somite_present_2 = get_somite_end(coord_2_trsf)

    somite_to_look = {}
    for t in set(somite_present_1.keys() + somite_present_2.keys()):
        somite_to_look[t] = min(max(somite_present_1.get(t, [somite_to_look.get(t-1, 2)])),
                                max(somite_present_2.get(t, [somite_to_look.get(t-1, 2)])))
    
    end_s_interp_1 = get_somite_interp(coord_1, somite_to_look)
    end_s_interp_2 = get_somite_interp(coord_2_trsf, somite_to_look)

    return end_s_interp_1, end_s_interp_2


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    if (v1 == 0).all() or (v2 == 0).all():
        return 0
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def oriented_angle_between(v1, v2, n):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    if (v1 == 0).all() or (v2 == 0).all():
        return 0
    else:
        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        A = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        if np.dot(n, np.cross(v1_u, v2_u)) > 0:
            return -A
        else:
            return A
        
def get_angles(params):
    ''' Calculates the angles between the equivalent somites
        between the two embryos
    '''
    t, origin = params
    try:
        start_time = time()
        b, n = symmetry_plane
        angles = []
        to_consider = []
        proj_ori = projection_onto_P(origin, b, n)[0] - b
        vE_2 = projection_onto_P(bary, b, n)[0] - b
        vB_2 = projection_onto_P(bary, b, n)[0] - b
        for i, L_2 in enumerate(LM_avg_rot_used):
            vP_2 = projection_onto_P(L_2, b, n)[0] - b
            a = oriented_angle_between(proj_ori, vP_2, n)
            if a < 0:
                a = 2 * np.pi - np.abs(a)
            if (np.linalg.norm(vP_2 - vE_2)>1 and 
                (angles==[] or (angles!=[] and np.round(a, 3) != angles[-1]))):
                to_consider += [i]
            angles += [np.round(a, 3)]
        angles = np.array(angles)[to_consider]
        tmp_LM_2 = [np.array(LM_avg_rot_used)[to_consider][i] for i in np.argsort(angles)]
        tmp_LM_1 = [np.array(LM_ref_rot)[to_consider][i] for i in np.argsort(angles)]
        angles = np.sort(angles)
        theta_disc = []
        for L_1, L_2 in zip(tmp_LM_1, tmp_LM_2):
            vP_1 = projection_onto_P(L_1, b, n)[0] - b
            vP_2 = projection_onto_P(L_2, b, n)[0] - b
            theta_disc += [oriented_angle_between(vP_1, vP_2, n)]
        theta_disc = np.array(theta_disc)
        return t, angles, theta_disc, to_consider, proj_ori
    except Exception as e:
        print e, t
        return(e, t)

def apply_rot(params):
    ''' Applies the piecewise rotation to the floating embryo
    '''
    t, angles, theta_disc, proj_ori, to_consider = params
    try:
        new_pos = {}
        theta = sp.interpolate.InterpolatedUnivariateSpline(angles,
                                                            theta_disc, k=1, ext = 3)
        b, n = symmetry_plane
        start_time = time()
        new_embryo = []
        previous_embryo = []
        new_embryo_angle = {}
        old_embryo_angle = {}
        for c in LT_for_multiprocess:
            pos = np.array(LT_for_multiprocess[c])
            previous_embryo += [pos]
            proj_pos = projection_onto_P(pos, b, n)[0] - b
            A_to_ORI = oriented_angle_between(proj_ori, proj_pos, n)
            if A_to_ORI < 0:
                A_to_ORI = 2 * np.pi - np.abs(A_to_ORI)
            theta_p = theta(A_to_ORI)
            T_v = rotation(pos - b, n, theta_p) + b
            new_pos[c] = T_v
            new_embryo += [T_v]
            new_embryo_angle[tuple(T_v)] = theta_p
            old_embryo_angle[tuple(pos)] = A_to_ORI

        tmp_LM_2 = [np.array(LM_avg_rot_used)[to_consider][i] for i in np.argsort(angles)]
        tmp_LM_1 = [np.array(LM_ref_rot)[to_consider][i] for i in np.argsort(angles)]
        LM_TRSF = []
        for L in tmp_LM_2:
            pos = L
            proj_pos = projection_onto_P(pos, b, n)[0] - b
            A_to_ORI = oriented_angle_between(proj_ori, proj_pos, n)
            if A_to_ORI < 0:
                A_to_ORI = 2 * np.pi - np.abs(A_to_ORI)
            theta_p = theta(A_to_ORI)
            T_v = rotation(pos - b, n, theta_p) + b
            LM_TRSF += [T_v]
        print 't %03d done in %.3f seconds'%(t, time() - start_time)
        return new_pos
    except Exception as e:
        print e, t
        return(e, t)


def apply_rot_P(params):
    ''' Applies the piecewise rotation to the floating embryo
    '''
    t, angles, theta_disc, proj_ori, to_consider, positions_given = params
    try:
        new_pos = {}
        theta = sp.interpolate.InterpolatedUnivariateSpline(angles,
                                                            theta_disc, k=1, ext = 3)
        b, n = symmetry_plane
        start_time = time()
        new_embryo = []
        previous_embryo = []
        new_embryo_angle = {}
        old_embryo_angle = {}
        for c in positions_given:
            pos = np.array(positions_given[c])
            previous_embryo += [pos]
            proj_pos = projection_onto_P(pos, b, n)[0] - b
            A_to_ORI = oriented_angle_between(proj_ori, proj_pos, n)
            if A_to_ORI < 0:
                A_to_ORI = 2 * np.pi - np.abs(A_to_ORI)
            theta_p = theta(A_to_ORI)
            T_v = rotation(pos - b, n, theta_p) + b
            new_pos[c] = T_v
            new_embryo += [T_v]
            new_embryo_angle[tuple(T_v)] = theta_p
            old_embryo_angle[tuple(pos)] = A_to_ORI

        tmp_LM_2 = [np.array(LM_avg_rot_used)[to_consider][i] for i in np.argsort(angles)]
        tmp_LM_1 = [np.array(LM_ref_rot)[to_consider][i] for i in np.argsort(angles)]
        LM_TRSF = []
        for L in tmp_LM_2:
            pos = L
            proj_pos = projection_onto_P(pos, b, n)[0] - b
            A_to_ORI = oriented_angle_between(proj_ori, proj_pos, n)
            if A_to_ORI < 0:
                A_to_ORI = 2 * np.pi - np.abs(A_to_ORI)
            theta_p = theta(A_to_ORI)
            T_v = rotation(pos - b, n, theta_p) + b
            LM_TRSF += [T_v]
        print 't %03d done in %.3f seconds'%(t, time() - start_time)
        return new_pos
    except Exception as e:
        print e, t
        return(e, t)


def dict_cmp(d1, d2):
    ''' Function to order two dictionaries
    '''
    return int(np.round(min(d1[0]) - min(d2[0])))

def extrapolate_pos(all_dicts, decay_time = 10):
    ''' Inter/extrapolate in time the position of the landmarks
        Args:
            all_dicts: {string: [float, float, float], }, dictionary that maps the name
                of a landmark to its position
    '''
    for tmp in sorted(all_dicts, cmp = dict_cmp):
        D, F = tmp
        t_start = min(D)
        if min(LT_2.new_time_nodes) < t_start:
            closest = np.inf
            for di, fi in all_dicts:
                times = [t for t in di.keys() if t <= t_start - decay_time]
                if times != []:
                    closest_t = max(times)
                    if np.linalg.norm(di[closest_t] - D[t_start]) < closest:
                        closest_dict = di
                        closest_f = fi
            D[t_start - decay_time] = closest_f(t_start - decay_time)
            for t, p in closest_dict.iteritems():
                if t < t_start - decay_time:
                    D[t] = p
        tmp[-1] = interp_3d_coord(np.array([np.array(list(v)+[k]) for k, v in D.iteritems()]))

    final_dists = []
    for D, F in all_dicts:
        final_dists += [interp_3d_coord(np.array([np.array(list(v)+[k]) for k, v in D.iteritems()]))]
    return final_dists

def compute_plane(p1, p2, p3):
    ''' Compute the parameters a, b, c, d, of a plane equation:
        ax + by + cz + d = 0 from 3 non-aligned points
        Args:
            p1, p2, p3: [float, float, float], 3D coordinates
        Returns:
            a, b, c, d: float, parameters of the equation of the plan
    '''
    v1 = p1 - p2
    v2 = p1 - p3
    cross_prod = np.cross(v1, v2)
    a, b, c = cross_prod
    d = np.sum([a, b, c] * -p1)
    return np.array([a, b, c, d])

def get_barycenter(fname):
    ''' Reads and coes a linear piecewise interpolation/extrapolation barycenters
        Args:
            fname: string, name of the barycenter file (each line as 'x, y, z, t')
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
    for l in lines[1:]:
        split_l = l.split(',')
        barycenters[int(split_l[-1])] = tuple(float(v) for v in split_l[:-1])
    times = sorted(barycenters)
    Xb, Yb, Zb = np.array([barycenters[t] for t in times]).T
    Xb_f = interpolate.InterpolatedUnivariateSpline(times, Xb, k=1)
    Yb_f = interpolate.InterpolatedUnivariateSpline(times, Yb, k=1)
    Zb_f = interpolate.InterpolatedUnivariateSpline(times, Zb, k=1)
    Ti = np.arange(11, 522)
    barycenters_interp = dict(zip(Ti, zip(Xb_f(Ti), Yb_f(Ti), Zb_f(Ti))))

    return barycenters_interp, barycenters

def build_gg(data):
    ''' Build a Gabriel graph from a list of 3D points.
        In the resulting Gabriel graph, the id of the nodes correspond
        to the position of points in *data*
        Args:
            data: Nx3 array of floats, array of 3D coordinates
        Returns:
            data: Nx3 array of floats, array of 3D coordinates (same as input)
            GG_out: {int: [int, ], }, dictionary that maps object id to the 
                list of the ids of its neighbors
            idx3d: kdtree structure
    '''
    D_graph = Delaunay(data, incremental = True)
    idx3d = kdtree(data)
    delaunay_graph = {}
    for N in D_graph.simplices:
        for e1, e2 in combinations(np.sort(N), 2):
            delaunay_graph.setdefault(e1, set([])).add(e2)
            delaunay_graph.setdefault(e2, set([])).add(e1)

    tmp = time()
    GG_out = {}
    for e1, neighbs in delaunay_graph.iteritems():
        for ni in neighbs:
            if not any([np.linalg.norm((data[ni] + data[e1])/2 - data[i])<np.linalg.norm(data[ni] - data[e1])/2
                    for i in delaunay_graph[e1].intersection(delaunay_graph[ni])]):
                GG_out.setdefault(e1, set()).add(ni)
                GG_out.setdefault(ni, set()).add(e1)
    return data, GG_out, idx3d

def get_V_W(P, GG, data, points1, pos):
    ''' From a set of 3D points, computes the set of vector+weights associated to
        including the vectors of the neighbors in the Gabriel graph
        Args:
            P; Nx3 array of floats, list of 3D points
            GG: {int: [int, ], }, dictionary that maps object id to the 
                list of the ids of its neighbors
            data: Mx3 array of floats, array of 3D points
            points1: Kx3 array of floats, array of 3D points
            pos: [float, float, float], 3D coordinate
        Returns:
            vector: Lx3 array (as a list), list of 3D vectors
            weights: Lx1 array (as a list), list of weights
    '''
    N = set()
    for pi in P:
        N.add(pi)
        for k in GG.get(pi, []):
            N.add(k)
    N = list(N)
    if N != []:
        weights = 1./np.linalg.norm(data[N] - pos, axis = 1)
        vector = points1[N] - data[N]
    return list(vector), list(weights)

def apply_trsf(t, k = 2):
    ''' Computes and apply a non-linear transformation between two embryos
        provided that the embryo shells Have been already computed.
        Args:
            t: int, time
            k: int, number of neighbors to consider (default 2)
    '''
    tic = time()
    points2_NL_1 = np.loadtxt(match_points_folder + 'pts02_low_DS_%03d.txt'%(t)) * DS_value
    points2_NL_2 = np.loadtxt(match_points_folder + 'pts02_high_DS_%03d.txt'%(t)) * DS_value

    points1_NL_1 = np.loadtxt(match_points_folder + 'pts01_low_DS_%03d.txt'%(t)) * DS_value
    points1_NL_2 = np.loadtxt(match_points_folder + 'pts01_high_DS_%03d.txt'%(t)) * DS_value

    data_1, GG_1, idx3d_1 = build_gg(points2_NL_1)
    data_2, GG_2, idx3d_2 = build_gg(points2_NL_2)
    neighbs = {}
    count = 0
    cells = new_pos_avg_SVF_rot.keys()
    closest_1 = idx3d_1.query([new_pos_avg_SVF_rot[c] for c in cells], k)
    closest_2 = idx3d_2.query([new_pos_avg_SVF_rot[c] for c in cells], k)
    final_pos = {}

    for i, (D, P) in enumerate(zip(*closest_1)):
        pos = new_pos_avg_SVF_rot[cells[i]]
        V_1, W_1 = get_V_W(P, GG_1, data_1, points1_NL_1, pos)
        V_2, W_2 = get_V_W(closest_2[1][i], GG_2, data_2, points1_NL_2, pos)

        V = np.array(V_1 + V_2)
        W = W_1 + W_2

        vector = np.sum(V * zip(W, W, W), axis = 0) / np.sum(W)

        final_pos[cells[i]] = pos + vector

    print 't%03d: %.2f'%(t, time() - tic)
    return final_pos

def apply_trsf_P(positions_given, t, k = 2):
    ''' Computes and apply a non-linear transformation between two embryos
        provided that the embryo shells Have been already computed.
        Args:
            t: int, time
            k: int, number of neighbors to consider (default 2)
    '''
    tic = time()
    points2_NL_1 = np.loadtxt(match_points_folder + 'pts02_low_DS_%03d.txt'%(t)) * DS_value
    points2_NL_2 = np.loadtxt(match_points_folder + 'pts02_high_DS_%03d.txt'%(t)) * DS_value

    points1_NL_1 = np.loadtxt(match_points_folder + 'pts01_low_DS_%03d.txt'%(t)) * DS_value
    points1_NL_2 = np.loadtxt(match_points_folder + 'pts01_high_DS_%03d.txt'%(t)) * DS_value

    data_1, GG_1, idx3d_1 = build_gg(points2_NL_1)
    data_2, GG_2, idx3d_2 = build_gg(points2_NL_2)
    neighbs = {}
    count = 0
    cells = positions_given.keys()
    closest_1 = idx3d_1.query([positions_given[c] for c in cells], k)
    closest_2 = idx3d_2.query([positions_given[c] for c in cells], k)
    final_pos = {}

    for i, (D, P) in enumerate(zip(*closest_1)):
        pos = positions_given[cells[i]]
        V_1, W_1 = get_V_W(P, GG_1, data_1, points1_NL_1, pos)
        V_2, W_2 = get_V_W(closest_2[1][i], GG_2, data_2, points1_NL_2, pos)

        V = np.array(V_1 + V_2)
        W = W_1 + W_2

        vector = np.sum(V * zip(W, W, W), axis = 0) / np.sum(W)

        final_pos[cells[i]] = pos + vector

    print 't%03d: %.2f'%(t, time() - tic)
    return final_pos

def rootname(n):
    new_n = n
    while new_n[-1].isdigit() or new_n[-1] == '.':
        new_n = new_n[:-1]
        if new_n[-1] == '.':
            new_n = new_n[:-1]
            break
    return new_n

def stabilize_point_cloud(SVF):
    SVF.compute_spatial_density(min(SVF.time_nodes), max(SVF.time_nodes))
    t_median = np.median(list(SVF.time_nodes))

    new_pos = {}
    for t in SVF.time_nodes:
        if t != t_median:
            if t < t_median:
                tree = SVF.successor
            elif t_median < t:
                tree = SVF.predecessor
            # else:
            #     break
            cells = SVF.time_nodes[t]
            tracks_to_ref = [] 
            for c in cells:
                final_c = c
                while SVF.time.get(final_c, t_median) != t_median:
                    final_c = tree.get(final_c, [-1])[0]
                if final_c != -1:
                    tracks_to_ref += [(c, final_c)]

            dens = np.array([SVF.spatial_density[c[0]] for c in tracks_to_ref])
            vect = np.array([SVF.pos[c[1]] - SVF.pos[c[0]] for c in tracks_to_ref])
            vect = np.apply_along_axis(lambda x:x*dens, 0, vect)
            avg_trsf = np.sum(vect, axis = 0)/np.sum(dens)
            new_pos.update({c: SVF.pos[c] + avg_trsf for c in cells})

    cells = SVF.time_nodes[t_median]
    new_pos.update({c:SVF.pos[c] for c in cells})
    return new_pos

def get_spherical_coordinates(x, y, z):
    ''' Computes spherical coordinates for an x, y, z Cartesian position
    '''
    r = np.linalg.norm([x, y, z])
    theta = np.arctan2(y, x)
    phi = np.arccos(z/r)
    alpha = (np.pi/2 + np.arctan2(x, z)) % (2*np.pi)
    return r, theta, phi, alpha

def write_DB(path_DB, path_div, VF, barycenters, tracking_value, tb, te):
    ''' Write the csv database in Database.csv
        Args:
            path_DB: string, path to the output database
            VF: lineageTree
            tracking_value: {int: int, }, dictionary that maps an object id to a label
            tb: int, first time point to write
            te: int, last time point to write
    '''
    f2 = open(path_DB + 'Database.csv', 'w')
    f2.write('id, mother_id, x, y, z, r, theta, phi, t, label, D-x, D-y, D-z, D-r, D-theta, D-phi\n')
    for t in range(tb, te+1):
        for c in VF.time_nodes[t]:
            S_p = (-1, -1, -1)
            if VF.predecessor.get(c, []) != []:
                M_id = VF.predecessor[c][0]
            else:
                M_id = -1
            P = tuple(VF.pos[c])
            S_p = tuple(get_spherical_coordinates(*(barycenters[t] - VF.pos[c]))[:-1])
            L = tracking_value.get(c, -1)
            D_P = (-1, -1, -1)
            D_S_p = (-1, -1, -1)
            f2.write(('%d, %d, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %d, %d,' + 
                      '%.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n')%((c, M_id) + P + S_p + (t, L) + D_P + D_S_p))
    f2.close()


def read_param_file():
    ''' Asks for, reads and formats the parameter file
    '''
    p_param = raw_input('Please enter the path to the parameter file:\n')
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
            if param_name in ['bary', 'DS_mask', 'labels']:
                name = param_name
                out = []
                while (name == param_name or param_name == '') and  i < nb_lines:
                    out += [float(split_line[1])]
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
        path_avg_TGMM = param_dict.get('path_avg_TGMM', '.')
        path_avg_SVF = param_dict.get('path_avg_SVF', '.')
        path_ref_TGMM = param_dict.get('path_ref_TGMM', '.')
        path_avg_LM = param_dict.get('path_avg_LM', '.')
        path_ref_LM = param_dict.get('path_ref_LM', '.')
        path_mask = param_dict.get('path_mask', '.')
        match_points_folder = param_dict.get('match_points_folder', '.')
        folder_SVF_amira = param_dict.get('folder_SVF_amira', '.')
        t_ref = int(param_dict.get('t_ref', 0))
        labels = np.array(param_dict.get('labels', [0])).astype(int)
        bary = np.array(param_dict.get('bary', [0, 0, 0]))
        transpose = int(param_dict.get('transpose_mask', 0))!=0
        DS_mask = np.array(param_dict.get('DS_mask', [1, 1, 1]))

    return (path_avg_TGMM, path_avg_SVF, path_ref_TGMM, path_avg_LM, path_ref_LM,
            match_points_folder, folder_SVF_amira, t_ref, labels, bary, path_mask, transpose, DS_mask)

if __name__ == '__main__':
    (path_avg_TGMM, path_avg_SVF, path_ref_TGMM, path_avg_LM,
        path_ref_LM, match_points_folder, folder_SVF_amira,
        t_ref, labels, bary, path_mask, transpose, DS_mask) = read_param_file()

    ### Preparation of the necessary folders
    ### TO DO ###
    if not os.path.exists(folder_SVF_amira):
        os.makedirs(folder_SVF_amira)
    if not os.path.exists(match_points_folder):
        os.makedirs(match_points_folder)

    corres_LM = {'HFAnt': 'HFAP',
                 'HFPost': 'HFPP',
                 'AntNot': 'AMN',
                 'SomL': 'LS',
                 'SomR': 'RS'}

    ### Reading of the different lineage trees
    # print path_avg_TGMM == '/data/Mouse-Project/algorithms/Average-embryo/initial_pos.bin'
    # print path_avg_SVF == '/data/Mouse-Project/algorithms/Average-embryo/avg_SVF_final.bin'
    # print path_ref_TGMM == '/media/R/SV1/KM_16-06-16/TGMM/GMEMtracking3D_ECC9_CV3.1/XML_finalResult_lht/GMEMfinalResult_frame%04d.xml'
    # print path_avg_LM == '/data/Mouse-Project/algorithms/Average-embryo/Landmarks/landmarks-avg-new.pkl'
    # print path_ref_LM == '/data/Mouse-Project/algorithms/Average-embryo/Landmarks/160616_shifted_TP265_ManualAnnotations.2.0.xml'
    # print match_points_folder == '/data/Mouse-Project/algorithms/Average-embryo/matching_points/'
    # print folder_SVF_amira == '/data/Shared_folder/T_amira/'
    # print t_ref == 265
    # print label == 250
    # print bary == np.array([ 1227.05137638,  1269.16905011,  1480.87914669])
    # bary[1] -= 200

    avg_TGMM = lineageTree(path_avg_TGMM)
    avg_SVF = lineageTree(path_avg_SVF)
    ref_TGMM = lineageTree(path_ref_TGMM, tb = t_ref, te = t_ref, z_mult = 5.)

    ### Reading of the landmarks
    with open(path_avg_LM) as f:
        LM_avg_rigid, LM_avg_rot = pkl.load(f)
    LM_avg_rigid = {em:{t:lms for t, lms in v.iteritems() if 2<len(lms)} for em, v in LM_avg_rigid.iteritems()}


    LM_ref = lineageTree(path_ref_LM, MaMuT = True)
    LM_ref.node_name = {k:corres_LM.get(v, v).replace('SOML', 'LS').replace('SOMR', 'RS')
                                                 for k, v in LM_ref.node_name.iteritems()}

    ### Finding the corresponding time (ie same number of somites)
    nb_som = len([v for v in LM_ref.node_name.values() if 'LS' == v[:2] or 'RS' == v[:2]])/2
    corres_time = np.round(np.mean([min(LM_avg_rot['RS%d'%nb_som]), min(LM_avg_rot['LS%d'%nb_som])]))

    ### Builds the first rigid transformation
    p2 = []
    p1 = []
    for LM_name, t_pos in LM_avg_rot.iteritems():
        if LM_name in LM_ref.node_name.values() and LM_name != 'ALM':
            keys = np.array(t_pos.keys())
            p2 += [t_pos[keys[np.argmin(np.abs(keys - corres_time))]]]
            p1 += [np.mean([LM_ref.pos[c] for c, n in LM_ref.node_name.iteritems() if n == LM_name], axis = 0)]

    mat = rigid_transform_3D(np.array(p2), np.array(p1))
    applyMatrixTrsf = lambda x, mat: list(np.dot(mat, np.transpose([list(x) + [1]])).reshape(4)[:3])

    ### Builds the second rigid transformation
    residual_best = np.inf
    for i, LMs in LM_avg_rigid.iteritems():
        keys = np.array(LMs.keys())
        LM_considered = LMs[keys[np.argmin(np.abs(keys - corres_time))]]
        final_mapping = {}
        for LM_name, positions in LM_considered.iteritems():
            corres_pos = [LM_ref.pos[c] for c, n in LM_ref.node_name.iteritems() if rootname(n) == LM_name]
            if corres_pos != []:
                final_mapping[LM_name] = [corres_pos, positions]
        points1 = []
        points2 = []
        for annotation_name, (cells1, cells2) in final_mapping.iteritems():
            a_ = np.array(cells1)
            b_ = np.array([applyMatrixTrsf(p, mat) for p in cells2])
            b_post_trsf = np.array(cells2)
            cost_matrix = spatial.distance_matrix(a_, b_)
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            points1 += list(a_[row_ind])
            points2 += list(b_post_trsf[col_ind])
        rigid_trsf = rigid_transform_3D(points2, points1)
        residual = np.mean([np.linalg.norm(applyMatrixTrsf(p2, rigid_trsf) - p1)
                                for p2, p1 in zip(points2, points1)])
        if residual < residual_best and 4 < len(LM_considered):
            trsf_best = rigid_trsf
            residual_best = residual
            i_final = i
            nb_LM = (len(LM_considered), np.sum([len(v) for v in LM_considered.iteritems()]))


    new_pos_avg_SVF = {c:applyMatrixTrsf(avg_SVF.pos[c], rigid_trsf) for c in avg_SVF.time_nodes[corres_time]}
    new_pos_avg_TGMM = {c:applyMatrixTrsf(avg_TGMM.pos[c], rigid_trsf) for c in avg_TGMM.time_nodes[corres_time]}

    DS_value = 10.
    VF_dimension = np.max(new_pos_avg_TGMM.values() + ref_TGMM.pos.values(), axis = 0)/DS_value + 50
    VF_dimension = np.round(VF_dimension)

    coord_ref = {LM_ref.node_name[c]: np.array(list(LM_ref.pos[c]) + [t_ref]) for c in LM_ref.time_nodes[t_ref]}
    line, somites_pairing = get_somite_pairing(coord_ref)

    avg = np.mean(line, axis = 0)
    uu, dd, vv = np.linalg.svd(line - avg)
    fitted_line = avg + vv[0] * (line[:, 0].reshape(-1, 1) - avg[0])

    projection_onto_P = lambda p, origin, n: ((p - np.dot(n, p - origin)/(np.linalg.norm(n))
                                              * (n/np.linalg.norm(n))), np.dot(n, p - origin)/(np.linalg.norm(n)))

    rotation = lambda v, r, theta: (1 - np.cos(theta)) * (np.dot(v, r)) * r + np.cos(theta) * v + np.sin(theta) * (np.cross(r, v))


    plane_params = compute_plane(np.array(bary), fitted_line[0], fitted_line[-1])
    symmetry_plane = bary, plane_params[:-1]/np.linalg.norm(plane_params[:-1])
    symmetry_plane_eq = plane_params
    coord_ref['PPT'] = coord_ref.get('PPT2', coord_ref['PPT'])

    LM_ref_rot = ([bary, coord_ref['HFPP'][:-1], coord_ref['AMN'][:-1]] +
                     [np.mean([coord_ref['RS%d'%i][:-1], coord_ref['LS%d'%i][:-1]], axis = 0)
                                        for i in range(1, nb_som+1)] +
                    [coord_ref['PPT2'][:-1]])

    LM_ref_rot_to_compare = ([coord_ref['HFPP'][:-1], coord_ref['AMN'][:-1]] + 
                             sum([[coord_ref['RS%d'%i][:-1], coord_ref['LS%d'%i][:-1]]
                                    for i in range(1, nb_som+1)], []) +
                             [coord_ref['PPT2'][:-1]])

    LM_avg_rot_used = [bary, ]
    LM_avg_rot_to_trsf = []
    for LM_n in [['HFPP'], ['AMN']] + [['RS%d'%i, 'LS%d'%i] for i in range(1, nb_som+1)] + [['PPT']]:
        pos_lm = []
        for lm_n in LM_n:
            keys = np.array(LM_avg_rot[lm_n].keys())
            pos_lm += [applyMatrixTrsf(LM_avg_rot[lm_n][keys[np.argmin(np.abs(keys - corres_time))]], rigid_trsf)]
        LM_avg_rot_used += [np.mean(pos_lm, axis = 0)]
        LM_avg_rot_to_trsf += pos_lm
        # LM_avg_rot_used += [pos_lm]
    # LM_ref_rot = ([bary, coord_ref['HFPP'][:-1], coord_ref['AMN'][:-1]] +
    #                  [np.mean([coord_ref['RS%d'%i][:-1], coord_ref['LS%d'%i][:-1]], axis = 0)
    #                                     for i in range(1, nb_som+1)])

    # LM_avg_rot_used = [bary, ]
    # for LM_n in [['HFPP'], ['AMN']] + [['RS%d'%i, 'RS%d'%i] for i in range(1, nb_som+1)]:
    #     pos_lm = []
    #     for lm_n in LM_n:
    #         keys = np.array(LM_avg_rot[lm_n].keys())
    #         pos_lm += [applyMatrixTrsf(LM_avg_rot[lm_n][keys[np.argmin(np.abs(keys - corres_time))]], rigid_trsf)]
    #     LM_avg_rot_used += [np.mean(pos_lm, axis = 0)]
    #     # LM_avg_rot_used += [pos_lm]

    origin = np.array([ 4000, 3000, 2500])
    t, angles, theta_disc, to_consider, proj_ori = get_angles([t_ref, origin])

    LT_for_multiprocess = new_pos_avg_TGMM
    new_pos_avg_TGMM_rot = apply_rot([t_ref, angles, theta_disc, proj_ori, to_consider])

    LT_for_multiprocess = new_pos_avg_SVF
    new_pos_avg_SVF_rot = apply_rot([t_ref, angles, theta_disc, proj_ori, to_consider])


    x, y, z = 0, 2, 1
    B_all = []
    R = 500
    nb_std = 1
    percentile = 80
    build_points(t_ref, bary)

    ## Computation of the distances
    # LM_avg_pos_dict = dict(zip(range(len(LM_avg_rot_to_trsf)), LM_avg_rot_to_trsf))
    # new_P = apply_rot_P([t_ref, angles, theta_disc, proj_ori, to_consider, LM_avg_pos_dict])
    # F_P = apply_trsf_P(new_P, t_ref)
    # distances = [np.linalg.norm(F_P[i] - LM_ref_rot_to_compare[i]) for i in range(len(F_P))]


    final_pos = apply_trsf(t_ref)

    ### Reading of the mask
    mask_dir = '.mask_images_tmp/'
    if not os.path.exists(mask_dir):
        os.makedirs(mask_dir)
    im = imread(path_mask)
    for l in labels:
        masked_im = im == l
        tmp = nd.binary_opening(masked_im, iterations = 3)
        tmp = nd.binary_closing(tmp, iterations = 4)
        imsave(mask_dir + '%03d.tif'%l, SpatialImage(tmp).astype(np.uint8))

    masks = [mask_dir + '%03d.tif'%l for l in labels]

    init_cells = {m: set() for m in range(len(masks))}

    for i, path_mask in enumerate(masks):
        if transpose:
            mask = imread(path_mask).transpose(1, 0, 2)
            mask = mask[:,::-1,:]
        else:
            mask = imread(path_mask)
        im_shape = mask.shape
        for c, p in final_pos.iteritems():
            pos_in_vox = np.round(p/DS_mask).astype(np.int)
            pos_in_vox = tuple(max(min(v, im_shape[i]-1), 0) for i, v in enumerate(pos_in_vox))
            if mask[pos_in_vox]:
                init_cells[i].add(c)

    tracking_value = {}
    for t, cs in init_cells.iteritems():
        for c in cs:
            to_treat = [c]
            tracking_value.setdefault(c, set()).add(t)
            while to_treat != []:
                c_tmp = to_treat.pop()
                next_cells = avg_SVF.successor.get(c_tmp, [])
                to_treat += next_cells
                for n in next_cells:
                    tracking_value.setdefault(n, set()).add(t)
            to_treat = [c]
            tracking_value.setdefault(c, set()).add(t)
            while to_treat != []:
                c_tmp = to_treat.pop()
                next_cells = avg_SVF.predecessor.get(c_tmp, [])
                to_treat += next_cells
                for n in next_cells:
                    tracking_value.setdefault(n, set()).add(t)

    tracking_value = {k:np.sum(list(v)) for k, v in tracking_value.iteritems() if len(v) == 1}

    # ref_mask = imread(path_mask)
    # masked_im = ref_mask == label
    # tmp = nd.binary_opening(masked_im, iterations = 3)
    # if transpose:
    #     mask = nd.binary_closing(tmp, iterations = 4)[::-1,:,:].transpose(1, 0, 2)
    # else:
    #     mask = nd.binary_closing(tmp, iterations = 4)
    # # mask[:,:,450:] = False

    # positive = set()
    # positive_pos = set()
    # im_shape = mask.shape
    # for c, p in final_pos.iteritems():
    #     pos_in_vox = tuple(np.round(p/DS_mask).astype(np.int))
    #     pos_in_vox = tuple(max(min(v, im_shape[i]-1), 0) for i, v in enumerate(pos_in_vox))
    #     if mask[pos_in_vox]:
    #         positive.add(c)
    #         positive_pos.add(tuple(p))

    # positive_dict = {}
    # for c in positive:
    #     track = [c]
    #     while track[-1] in avg_SVF.successor:
    #         track += avg_SVF.successor[track[-1]]
    #     while track[0] in avg_SVF.predecessor:
    #         track.insert(0, avg_SVF.predecessor[track[0]][0])
    #     positive_dict.update({k:1 for k in track})

    new_pos = stabilize_point_cloud(avg_SVF)

    # timeTrsf = sp.interpolate.InterpolatedUnivariateSpline([91, 252], [20, 265], k=1)
    # avg_SVF.pos = new_pos
    # tmp = apply_time_trsf(avg_SVF, timeTrsf, min_t = min(avg_SVF.time_nodes), max_t = max(avg_SVF.time_nodes))

    if not os.path.exists(folder_SVF_amira):
        os.makedirs(folder_SVF_amira)
    write_to_am_2(folder_SVF_amira + '/seg_t%04d.am', avg_SVF, t_b = None, t_e = None,
                    manual_labels = tracking_value, default_label = 0, length = 7,
                      new_pos = new_pos)#, to_take_time = {k+29:v for k, v in avg_SVF.new_time_nodes.iteritems()})

    for im_p in os.listdir(mask_dir):
        os.remove(mask_dir + '/' + im_p)

