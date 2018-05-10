# This file is subject to the terms and conditions defined in
# file 'LICENSE', which is part of this source code package.
# Author: Leo Guignard (guignardl...@AT@...janelia.hhmi.org)

import os
import sys
from time import time
import struct
from multiprocessing import Pool
from itertools import combinations
import xml.etree.ElementTree as ET
import numpy as np
from scipy.spatial import Delaunay
from scipy import spatial
import scipy as sp
from scipy import ndimage

from TGMMlibraries import lineageTree

def single_cell_propagation(params):
    ''' Computation of the median vector of a cell *C*.
        This function is suited for parallel computation.
    '''
    C, closest_cells, dist_max, W = params
    dists = np.array([np.sum((LT.pos[C]-LT.pos[n])**2)**.5 for n in closest_cells])
    if (dists<dist_max).any():
        cells_to_keep = closest_cells[np.where(dists<dist_max)]
        W_to_keep = W[np.where(dists<dist_max)]
        subset_dist = [np.mean([LT.pos[cii] for cii in LT.predecessor[ci]], axis=0) - LT.pos[ci] 
                       for ci in cells_to_keep if not LT.is_root[ci]]
        subset_W = [W_to_keep[i] for i, ci in enumerate(cells_to_keep) if not LT.is_root[ci]]
        if subset_dist != []:
            med_distance = spatial.distance.squareform(spatial.distance.pdist(subset_dist)) * subset_W
            med = subset_dist[np.argmin(np.sum(med_distance, axis=1))]
        else:
            med = np.array([0, 0, 0])
    else:
        med = np.array([0, 0, 0])
    return C, med

def build_VF_propagation_backward(LT, t_b=0, t_e=200, neighb_size=20, dist_max=200, nb_proc = 8, Ws = [4, 1, 0]):
    ''' Main function that performs the backward propagation.
        The SVF is saved in LT.VF
        Args:
            LT: lineageTree
            t_b: int, first time point to process (note that t_e < t_b since we are propagating backward)
            t_e: int, last time point to process
            neighb_size: int, number of neighbors to look at to build the neighborhood density
            dist_max: float, distance max from which a vector is considered.
                Notice that the value used is so big that in that case it does not discard any cell.
                This parameter is historical
            nb_proc: int, number of processes to run in parallel:
                1: no parallelisation
                negative number: as many as possible (default)
                1 < *nb*: *nb* parallel processes.
            Ws: [int, int, int], list of weights for the neighborhood
    '''
    if (not hasattr(LT, 'VF')) or LT.VF == None:
        LT.VF = lineageTree(None, None, None)
        starting_cells = LT.time_nodes[t_b]
        unique_id = 0
        LT.VF.time_nodes = {t_b: []}
        for C in starting_cells:
            i = LT.VF.get_next_id()
            LT.VF.nodes.append(i)
            LT.VF.time_nodes[t_b].append(i)
            LT.VF.roots.append(i)
            LT.VF.pos[i]=LT.pos[C]

    gg_line = '\r%03d: GG done (%.2f s); '
    prop_line = 'P done (%.2f s); '
    add_line = 'Add done (%.2f s); '
    fusion_line = 'F done (%.2f s); '
    nb_cells_line = '#C: %06d'

    full_length = len(gg_line) + 3 + 5 +\
                  len(prop_line) + 5 +\
                  len(add_line) + 5 + \
                  len(fusion_line) + 5 + \
                  len(nb_cells_line) + 6

    
    for t in range(t_b, t_e, -1):
        tic = time()

        to_check_VF = LT.VF.time_nodes[t]
        idx3d, to_check_LT = LT.get_idx3d(t)
        LT.VF.time_nodes[t-1] = []
        Gabriel_graph = LT.get_gabriel_graph(t)
        GG_pred = LT.get_gabriel_graph(t-1) if t_e < t else {}
        GG_succ = LT.get_gabriel_graph(t+1) if t < t_b else {}

        gg_time = time() - tic
        tagged = set()
        cell_mapping_LT_VF = {}
        mapping = []
        for C in to_check_VF:
            C_LT = to_check_LT[idx3d.query(LT.VF.pos[C])[1]]
            cell_mapping_LT_VF.setdefault(C_LT, []).append(C)
            if not C_LT in tagged:
                C_LT_pred = set(LT.predecessor.get(C_LT, []))
                C_LT_pred_2 = set([ci for pred in LT.predecessor.get(C_LT, []) for ci in GG_pred.get(pred, set())])
                C_LT_succ = set(LT.successor.get(C_LT, []))
                C_LT_succ_2 = set([ci for succ in LT.successor.get(C_LT, []) for ci in GG_succ.get(succ, set())])
                N = list(Gabriel_graph.get(C_LT, set()))
                weight_mat = np.ones_like(N)
                N_pred = {n: LT.predecessor.get(n, []) for n in N}
                N_succ = {n: LT.successor.get(n, []) for n in N}
                for i, n in enumerate(N):
                    for n_predii in N_pred[n]:
                        if C_LT_pred.intersection(GG_pred.get(n_predii, set())):
                            W = Ws[0]
                        elif C_LT_pred_2.intersection(GG_pred.get(n_predii, set())):
                            W = Ws[1]
                        else:
                            W = Ws[2]
                    for n_succii in N_succ[n]:
                        if C_LT_succ.intersection(GG_succ.get(n_succii, set())):
                            W += Ws[0]
                        elif C_LT_succ_2.intersection(GG_succ.get(n_succii, set())):
                            W += Ws[1]
                        else:
                            W += Ws[2]
                    weight_mat[i] = W
                if (weight_mat == 0).all():
                    weight_mat[:] = 1
                mapping += [(C_LT, np.array(list(N)), dist_max, weight_mat)]
                tagged.add(C_LT)

        out = []

        if nb_proc == 1:
            for params in mapping:
                out += [single_cell_propagation(params)]
        else:
            if nb_proc < 1:
                pool = Pool()
            else:
                pool = Pool(processes=nb_proc)
            out = pool.map(single_cell_propagation, mapping)
            pool.terminate()
            pool.close()
        for C_LT, med in out:
            for C_VF in cell_mapping_LT_VF[C_LT]:
                LT.VF.add_node(t-1, C_VF, LT.VF.pos[C_VF] + med, reverse = True)

        idx3d, to_check_LT = LT.get_idx3d(t-1)
        to_check_VF = LT.VF.time_nodes[t-1]

        sys.stdout.write('\b'*(full_length) + ' '*(full_length))
        sys.stdout.flush()
        sys.stdout.write(gg_line%(t, gg_time))
        sys.stdout.flush()

        sys.stdout.write(prop_line%(time() - tic))
        sys.stdout.flush()

        if not LT.spatial_density.has_key(to_check_LT[0]):
            LT.compute_spatial_density(t-1, t-1, neighb_size)

        idx3d, to_check_VF = LT.VF.get_idx3d(t-1)[:2]
        dist_to_VF, equivalence = idx3d.query([LT.pos[c] for c in to_check_LT], 1)
        tmp = np.array([dist_to_VF[i]/LT.spatial_density[c] for i, c in enumerate(to_check_LT)])
        to_add = [to_check_LT[i] for i in np.where(tmp>1.25)[0]]
        for C in to_add:
            LT.VF.add_node(t-1, None, LT.pos[C])

        sys.stdout.write(add_line%(time() - tic))
        sys.stdout.flush()
        sys.stdout.write(nb_cells_line%i)
        sys.stdout.flush()

    LT.VF.t_b = t_e
    LT.VF.t_e = t_b

    return LT.VF

def prune_tracks(VF, mapping_LT_to_VF):
    ''' Prunes the SVF lineage tree *VF* by removing the objects in *VF* that are
        not consistantly mapped to objects in the TGMM lineage tree.
        Args:
            VF: lineageTree
            mapping_LT_to_VF: {int: [int, ] }, dictionary that maps a object in TGMM
                    onto its 10 closest objects in the SVF
    '''
    to_keep = {}
    for v in mapping_LT_to_VF.itervalues():
        for c in v:
            to_keep[c] = True

    to_keep_tmp = {}
    for c in to_keep.iterkeys():
        m = c
        d = c
        for i in range(5):
            m = VF.successor.get(m, [m])[0]
            d = VF.predecessor.get(d, [d])[0]
            to_keep_tmp[m] = True
            to_keep_tmp[d] = True
        to_keep_tmp[c] = True

    to_keep_final = {}
    for c in to_keep_tmp.iterkeys():
        m = c
        d = c
        i = 0
        while (    to_keep_tmp.get(VF.successor.get(m, [m])[0], False)
               and to_keep_tmp.get(VF.predecessor.get(d, [d])[0], False)
               and i < 5):
            m = VF.successor.get(m, [m])[0]
            d = VF.predecessor.get(d, [d])[0]
            i += 1

        to_keep_final[c] = (i == 5)

    VF.to_take_time = {}
    for C in VF.nodes:
        if to_keep_final.get(C, False) and C in VF.time:
            VF.to_take_time.setdefault(VF.time[C], []).append(C)

def get_gabriel_graph_for_parallel(params):
    ''' Computes the Gabriel graph for the global lineage tree *LT*
        at time point t.
    '''
    t = params
    if not hasattr(LT, 'Gabriel_graph'):
        LT.Gabriel_graph = {}

    if not LT.Gabriel_graph.has_key(t):
        idx3d, nodes = LT.get_idx3d(t)

        data_corres = {}
        data = []
        for i, C in enumerate(nodes):
            data.append(LT.pos[C])
            data_corres[i] = C

        tmp = Delaunay(data)

        delaunay_graph = {}

        for N in tmp.simplices:
            for e1, e2 in combinations(np.sort(N), 2):
                delaunay_graph.setdefault(e1, set([])).add(e2)
                delaunay_graph.setdefault(e2, set([])).add(e1)


        Gabriel_graph = {}

        for e1, neighbs in delaunay_graph.iteritems():
            for ni in neighbs:
                if not any([np.linalg.norm((data[ni] + data[e1])/2 - data[i])<np.linalg.norm(data[ni] - data[e1])/2
                        for i in delaunay_graph[e1].intersection(delaunay_graph[ni])]):
                    Gabriel_graph.setdefault(data_corres[e1], set()).add(data_corres[ni])
                    Gabriel_graph.setdefault(data_corres[ni], set()).add(data_corres[e1])
    else:
        Gabriel_graph = LT.Gabriel_graph[t]

    return t, Gabriel_graph

def parallel_gabriel_graph_preprocess(LT, nb_proc = 24):
    ''' Computes the gabriel graphs for each time point of a lineage tree in parallel.
        The results is written in LT.Gabriel_graph
        Args:
            LT: lineageTree
            nb_proc: int, number of processes to run in parallel:
                1: no parallelisation
                negative number: as many as possible (default)
                1 < *nb*: *nb* parallel processes.
    '''
    mapping = []
    if not hasattr(LT, 'Gabriel_graph'):
        LT.Gabriel_graph = {}
    for t in xrange(LT.t_b, LT.t_e + 1):
        if not LT.Gabriel_graph.has_key(t):
            mapping += [(t)]
    if nb_proc == 1:
        out = []
        for params in mapping:
          out += [get_gabriel_graph_for_parallel(params)]
    else:
        if nb_proc < 1:
            pool = Pool(processes)
        else:
            pool = Pool(processes=nb_proc)
        out = pool.map(get_gabriel_graph_for_parallel, mapping)
        pool.terminate()
        pool.close()
    for t, G_g in out:
        LT.Gabriel_graph[t] = G_g

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
                if 0 < len(LT_to_print.predecessor.get(C_tmp, [C_tmp])):
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

def mapping_VF_to_LT(params):
    ''' Maps objects in the global TGMM lineage tree *LT* to the 10 closest
        objects in the global SVF, *VF*, for time point *t*
        Suited for parallel processing (hence the single argument)
        Args:
            params: (t,), time to process.
        returns:
            final_mapping: {id1: id2, }, dictionary where the keys are node ids in VF
                            and the values are the associated node ids in LT for time point t
    '''
    t, = params
    cells_VF = VF.time_nodes[t]
    cells_LT = LT.time_nodes[t]
    idx3d, mapping_VF = VF.get_idx3d(t)

    positions_LT = [LT.pos[C] for C in cells_LT]
    mapping = idx3d.query(positions_LT, 10)[1]
    final_mapping = {}
    for i, C_neighs in enumerate(mapping):
        final_mapping[cells_LT[i]] = [mapping_VF[ci] for ci in C_neighs]
    return final_mapping

def parallel_mapping(tb, te, nb_proc = -1):
    ''' Maps objects in the global TGMM lineage tree *LT* to their closest objects in the global SVF, *VF*.
        Args:
            tb: int, starting time point
            te: int, ending time point
            nb_proc: int, number of processes to run in parallel:
                1: no parallelisation
                negative number: as many as possible (default)
                1 < *nb*: *nb* parallel processes.
        Returns:
            mapping_VF_to_LT_out: {id1: id2, }, dictionary where the keys are node ids in VF
                                    and the values are the associated node ids in LT
    '''
    mapping = []
    for t in range(tb, te+1):
        mapping += [[t]]

    if nb_proc == 1:
        out = []
        for params in mapping:
            out += [mapping_VF_to_LT(params)]
    else:
        if nb_proc < 1:
            pool = Pool()
        else:
            pool = Pool(processes = nb_proc)
        out = pool.map(mapping_VF_to_LT, mapping)
        pool.terminate()
        pool.close()

    mapping_VF_to_LT_out = {}
    for tmp in out:
        mapping_VF_to_LT_out.update(tmp)

    return mapping_VF_to_LT_out

def GG_to_bin(gg, fname):
    ''' Write a Gabriel graph to a binary file.
        The file is composed of a list of nodes, a list of couples (times, #nodes for that time)
        and a header which is the size of each of the two previously mentioned lists.
        The list of couples (time, #nodes for that time), *time_list*, stores the number of nodes
            enumerated in *nodes_list* that belong to that time point. For example, if the list is
            [(0, 120), (1, 341), (10, 432)] it means that the 120 first node ids in *nodes_list* belong to
            TP 0, the next 341 belong to TP 1 and the next 432 belong to TP 10.
        The list of nodes, *nodes_list*, is a list of node ids.
            A negative number *nn* means that this is the source id (-*nn*), all the positive node
            ids following that number are nodes link to -*nn*. The next negative number is a new source id.
            For example, [-1, 2, 3, -2, 1, 3, -3, 1, 2], would discribe the complete digraph K3.
        The header state the length of *time_list*, and the length of *nodes_list*.
        Each value is stored as long long (0 -> 2^(8*8)-1).
        Args:
            gg: {source: [target1, ...]}, dictionary with node source id as keys and a list of target ids as values.
            fname: string, file name
    '''
    nodes_list = []
    time_list = []
    for t, gg_t in gg.iteritems():
        time_list += [t]
        len_ = 0
        for node, neighbors in gg_t.iteritems():
            to_add = [-(node + 1)] + [neighb + 1 for neighb in neighbors]
            nodes_list += to_add
            len_ += len(to_add)
        time_list += [len_]

    f = open(fname, 'wb')
    f.write(struct.pack('q', len(time_list)))
    f.write(struct.pack('q', len(nodes_list)))
    f.write(struct.pack('q'*len(time_list), *time_list))
    f.write(struct.pack('q'*len(nodes_list), *nodes_list))
    f.close()


def GG_from_bin(fname):
    ''' Reads a Gabriel graph from a binary file
        Args:
            fname: string, name of the binary file
    '''
    q_size = struct.calcsize('q')
    f = open(fname)
    len_time_list = struct.unpack('q', f.read(q_size))[0]
    len_nodes_list = struct.unpack('q', f.read(q_size))[0]
    time_list = list(np.reshape(struct.unpack('q'*len_time_list, 
                                              f.read(q_size*len_time_list)),
                                (len_time_list/2, 2)))
    nodes_list = list(struct.unpack('q'*len_nodes_list, f.read(q_size*len_nodes_list)))
    f.close()

    gg ={}
    pos = 0
    for t, len_ in time_list:
        gg[t] = {}
        for n in nodes_list[pos:pos + len_]:
            if n < 0:
                current_node = -n - 1
                gg[t][current_node] = set()
            else:
                gg[t][current_node].add(n - 1)
        pos += len_

    return gg

def read_param_file():
    ''' Asks for, reads and formats the csv parameter file.
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
            if param_name in []:
                name = param_name
                out = []
                while (name == param_name or param_name == '') and  i < nb_lines:
                    out += [int(split_line[1])]
                    i += 1
                    if i < nb_lines:
                        l = lines[i]
                        split_line = l.split(',')
                        param_name = split_line[0]
                param_dict[name] = np.array(out)
            else:
                param_dict[param_name] = split_line[1].strip()
                i += 1
            if param_name == 'anisotropy' or 'time' in param_name:
                if split_line[1].isdigit():
                    param_dict[param_name] = int(split_line[1])
                else:
                    param_dict[param_name] = float(split_line[1])
        path_to_xml = param_dict.get('path_to_xml', '.')
        path_out = param_dict.get('path_out', '.')
        anisotropy = param_dict.get('anisotropy', 1)
        start_time = param_dict.get('start_time', -1)
        end_time = param_dict.get('end_time', -1)
        p_TGMM = param_dict.get('path_bin_TGMM', None)
    return (path_to_xml, path_out, anisotropy, start_time, end_time, p_TGMM)

if __name__ == '__main__':
    global LT
    path_to_xml, path_out, anisotropy, tb, te, p_TGMM = read_param_file()
    
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    if not os.path.exists(path_out + 'Amira_SVF/'):
        os.makedirs(path_out + 'Amira_SVF/')
    if not os.path.exists(path_out + 'Amira_TGMM/'):
        os.makedirs(path_out + 'Amira_TGMM/')

    if p_TGMM is None or not os.path.exists(p_TGMM):
        files = [f for f in os.listdir(path_to_xml) if '.xml' in f]
        pos_time = len(os.path.commonprefix(files))
        times = [int(file.split('.')[0][pos_time:]) for file in files]
        if tb < 0:
            tb = min(times) + 10
        if te < 0:
            te = max(times) - 10
        LT_main = lineageTree(file_format = path_to_xml + '/GMEMfinalResult_frame%04d.xml',
                         tb = tb, te = te, z_mult = anisotropy)
        LT_main.to_binary(path_out + 'TGMM.bin')
    else:
        LT_main = lineageTree(file_format = p_TGMM)
        tb = LT_main.t_b
        te = LT_main.t_e
    LT = LT_main

    if not os.path.exists(path_out + 'GG.bin'):
        LT = LT_main
        tic = time()
        parallel_gabriel_graph_preprocess(LT_main)
        GG_to_bin(LT_main.Gabriel_graph, path_out + 'GG.bin')
        print 'Gabriel graph pre-processing:',  time() - tic
    else:
        LT_main.Gabriel_graph = GG_from_bin(path_out + 'GG.bin')

    if not os.path.exists(path_out + 'SVF.bin'):
        tic = time()
        VF = build_VF_propagation_backward(LT_main, t_b = te, t_e = tb, neighb_size = 10, nb_proc=24)
        print 'parallel processing:',  time() - tic

        mapping_LT_to_VF = parallel_mapping(tb, te)
        prune_tracks(VF, mapping_LT_to_VF)

        done = set()
        corresponding_track = {}
        smoothed_pos = {}
        num_track = 0
        for C in VF.nodes:
            if not C in done:
                track = [C]
                while track[-1] in VF.successor:
                    track.append(VF.successor[track[-1]][0])
                while track[0] in VF.predecessor:
                    track.insert(0, VF.predecessor[track[0]][0])
                pos_track = np.array([VF.pos[Ci] for Ci in track])
                X = sp.ndimage.filters.gaussian_filter1d(pos_track[:, 0], sigma = 5)
                Y = sp.ndimage.filters.gaussian_filter1d(pos_track[:, 1], sigma = 5)
                Z = sp.ndimage.filters.gaussian_filter1d(pos_track[:, 2], sigma = 5)
                track_smoothed = np.zeros_like(pos_track)
                track_smoothed[:, 0] = X
                track_smoothed[:, 1] = Y
                track_smoothed[:, 2] = Z
                smoothed_pos.update(zip(track, list(track_smoothed)))
                done.update(set(track))
        
        to_remove = set()
        for t, n in VF.to_take_time.iteritems():
            to_remove.update(set(VF.time_nodes[t]).difference(n))

        VF.nodes = set(VF.nodes)
        for c in to_remove:
            if VF.predecessor.get(c, []) != [] and c in VF.successor.get(VF.predecessor.get(c, [-1])[0], []):
                VF.successor[VF.predecessor[c][0]].remove(c)
            for ci in VF.successor.get(c, []):
                if c in VF.predecessor.get(ci, []):
                    VF.predecessor.get(ci, []).remove(c)
            VF.successor.pop(c, [])
            VF.successor.pop(c, [])
            VF.predecessor.pop(c, [])
            VF.nodes.remove(c)

        VF.pos = smoothed_pos

        VF.to_binary(path_out + 'SVF.bin')
    else:
        VF = lineageTree(path_out + 'SVF.bin')

    write_to_am_2(path_out + 'Amira_SVF/seg_t%04d.am', VF, t_b = None, t_e = None,
                  manual_labels = {}, default_label = 1, length = 7)

    f2 = open(path_out + 'Database-SVF.csv', 'w')
    f2.write('id, mother_id, x, y, z, r, theta, phi, t, label, D-x, D-y, D-z, D-r, D-theta, D-phi\n')
    f1 = open(path_out + 'Database_expanded-SVF.csv', 'w')
    f1.write('id, mother_id, x, y, z, r, theta, phi, t, label, D-x, D-y, D-z, D-r, D-theta, D-phi,' +
             ' x_M, y_M, z_M, r_M, theta_m, phi_m, d_m_x, d_m_y, d_m_z, d_m_r, d_m_theta, d_m_phi,' +
             ' d_m_L2, alpha, M_alpha, d_m_alpha, D_alpha\n')
    path_bary = None
    tracking_value = {}
    ass_div = {}
    for t in range(tb, te+1):
        for c in VF.time_nodes[t]:
            S_p = (-1, -1, -1)
            alpha = -1
            S_p = (-1, -1, -1)
            m_S_p = (-1, -1, -1)
            m_alpha = -1
            d_m_alpha = -1
            d_m_P = (-1, -1, -1)
            d_m_S_p = (-1, -1, -1)
            d_m_L2 = -1
            m_P = (-1, -1, -1)
            if VF.predecessor.get(c, []) != []:
                M_id = VF.predecessor[c][0]
            else:
                M_id = -1
            P = tuple(VF.pos[c])
            if path_bary is not None:
                S_p = tuple(get_spherical_coordinates(*(barycenters[t] - VF.pos[c]))[:-1])
                alpha = get_spherical_coordinates(*(barycenters[t] - VF.pos[c]))[-1]
            if M_id != -1:
                m_P = tuple(VF.pos[M_id])
                d_m_P = tuple(np.array(P) - m_P)
                d_m_L2 = np.linalg.norm(d_m_P)
                if path_bary is not None:
                    m_S_p = tuple(get_spherical_coordinates(*(barycenters[t-1] - VF.pos[M_id]))[:-1])
                    m_alpha = get_spherical_coordinates(*(barycenters[t-1] - VF.pos[M_id]))[-1]
                    d_m_alpha = alpha - m_alpha
                    d_m_S_p = tuple(np.array(S_p) - m_S_p)
            L = tracking_value.get(c, -1)
            D_P = tuple(ass_div.get(c, [-1, -1, -1]))
            if path_bary is not None:
                D_S_p = (-1, -1, -1) if not c in ass_div else tuple(get_spherical_coordinates(*(barycenters[t] - ass_div[c]))[:-1])
                D_alpha = -1 if not c in ass_div else get_spherical_coordinates(*(barycenters[t] - ass_div[c]))[-1]
            else:
                D_S_p = (-1, -1, -1)
                D_alpha = -1
            f2.write(('%d, %d, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %d, %d,' + 
                      '%.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n')%((c, M_id) + P + S_p + (t, L) + D_P + D_S_p))
            f1.write(('%d, %d, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %d, %d,' + 
                      '%.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f,'+
                      '%.5f, %.5f, %.5f,%.5f, %.5f, %.5f, %.5f, %.5f, %.5f,' + 
                      '%.5f, %.5f, %.5f, %.5f, %.5f\n')%((c, M_id) + P + S_p + (t, L) + D_P + D_S_p + m_P +
                                                         m_S_p + d_m_P + d_m_S_p + (d_m_L2, alpha, m_alpha, d_m_alpha, D_alpha)))
    
    f1.close()
    f2.close()

    done = set()
    corresponding_track = {}
    smoothed_pos = {}
    num_track = 0
    all_tracks = []
    for C in LT.nodes:
        if not C in done:
            track = [C]
            while track[-1] in LT.successor:
                track.append(LT.successor[track[-1]][0])
            while track[0] in LT.predecessor:
                track.insert(0, LT.predecessor[track[0]][0])
            all_tracks += [track]
            done.update(set(track))
            pos_track = np.array([LT.pos[Ci] for Ci in track])
            X = sp.ndimage.filters.gaussian_filter1d(pos_track[:, 0], sigma = 5)
            Y = sp.ndimage.filters.gaussian_filter1d(pos_track[:, 1], sigma = 5)
            Z = sp.ndimage.filters.gaussian_filter1d(pos_track[:, 2], sigma = 5)
            track_smoothed = np.zeros_like(pos_track)
            track_smoothed[:, 0] = X
            track_smoothed[:, 1] = Y
            track_smoothed[:, 2] = Z
            smoothed_pos.update(zip(track, list(track_smoothed)))

    write_to_am_2(path_out + 'Amira_TGMM/seg_t%04d.am', LT, t_b = None, t_e = None,
                  manual_labels = {}, default_label = 1, 
                  length = 7, new_pos = smoothed_pos)

