from TGMMlibraries import lineageTree
import numpy as np
import os
from mamut_xml_templates import *


def read_param_file():
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
            if param_name in ['labels', 'downsampling', 'tissue_ids', 'tissue_names']:
                name = param_name
                out = []
                while (name == param_name or param_name == '') and  i < nb_lines:
                    if split_line[1].strip().isdigit():
                        out += [int(split_line[1])]
                    else:
                        out += [(split_line[1].strip())]
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
        path_SVF = param_dict.get('path_to_SVF', '.')
        path_DB = param_dict.get('path_to_DB', '.')
        path_output = param_dict.get('path_output', '.')
        tissue_ids = param_dict.get('tissue_ids', [])
        tissue_names = param_dict.get('tissue_names', [])
        begin = int(param_dict.get('begin', None))
        end = int(param_dict.get('end', None))
        filename = param_dict.get('filename', '.')
        folder = param_dict.get('folder', '.')
        v_size = np.float(param_dict.get('v_size', 1.))
        dT = np.float(param_dict.get('dT', 1.))

    return (path_SVF, path_DB, path_output, tissue_ids, 
        tissue_names, begin, end, filename, folder, v_size, dT)

def main():
    try:
        (path_SVF, path_to_DB, path_output, tissue_ids, 
            tissue_names, begin, end, filename, folder, v_size, dT) = read_param_file()
    except Exception as e:
        print "Failed at reading the configuration file."
        print "Error: %s"%e
        raise e
    
    SVF = lineageTree(path_SVF)

    # f = open(path_to_DB)
    # lines = f.readline()
    # f.close()

    DATA = np.loadtxt(path_to_DB, delimiter = ',', skiprows = 1, usecols = (0, 9))

    tracking_value = dict(DATA[:, (0, 1)])

    kept_times = range(begin, end+1)
    kept_nodes = [c for c in SVF.nodes if tracking_value[c] in tissue_ids]
    kept_nodes_set = set(kept_nodes)
    first_nodes = [c for c in SVF.time_nodes[min(kept_times)] if c in kept_nodes_set]
    t_id_2_N = dict(zip(tissue_ids, tissue_names))

    if not os.path.exists(os.path.dirname(path_output)):
        os.makedirs(os.path.dirname(path_output))
    with open(path_output, 'w') as output:

        output.write(begin_template)

        # Begin AllSpots.
        output.write(allspots_template.format(nspots=len(kept_nodes)))

        # Loop through lists of spots.
        for t in kept_times:
            cells = kept_nodes_set.intersection(SVF.time_nodes.get(t, []))
            if cells != []:
                output.write(inframe_template.format(frame=t))
                for c in cells:
                    output.write(spot_template.format(id=c, name=c, frame=t, t_id=tracking_value[c],
                                                      x=SVF.pos[c][0]*v_size,
                                                      y=SVF.pos[c][1]*v_size,
                                                      z=SVF.pos[c][2]*v_size,
                                                      # t_name=t_id_2_N.get(tracking_value[c], '?')
                                                      ))
                output.write(inframe_end_template)
            else:
                output.write(inframe_empty_template.format(frame=t))


        # End AllSpots.
        output.write(allspots_end_template)

        all_tracks = []
        roots = set(kept_nodes).difference(SVF.predecessor).union(first_nodes)
        last_time = max(kept_times)
        for c in roots:
            track = [c]
            while track[-1] in SVF.successor and SVF.time[track[-1]]<last_time:
                track += SVF.successor[track[-1]]
            all_tracks += [track]

        # Begin AllTracks.
        output.write(alltracks_template)

        for track_id, track in enumerate(all_tracks[::-1]):
            stop = SVF.time[track[-1]]
            duration = len(track)
            output.write(track_template.format(id=track_id+1, duration=duration, 
                                               start=SVF.time[track[0]], stop=stop, nspots=len(track),
                                               displacement=np.linalg.norm(SVF.pos[track[0]]-SVF.pos[track[-1]])))
            for c in track[:-1]:
                displacement = np.linalg.norm(SVF.pos[c] - SVF.pos[SVF.successor[c][0]]) * v_size
                velocity = displacement / dT
                output.write(edge_template.format(source_id=c, target_id=SVF.successor[c][0],
                                                  # t_name=t_id_2_N.get(tracking_value[c], '?'),
                                                  velocity=velocity, displacement=displacement,
                                                  t_id=tracking_value[c], time=SVF.time[c]))
            output.write(track_end_template)

        # End AllTracks.
        output.write(alltracks_end_template)

        # Filtered tracks.
        output.write(filteredtracks_start_template)
        for track_id, track in enumerate(all_tracks[::-1]):
            output.write(filteredtracks_template.format(t_id=track_id+1))
        output.write(filteredtracks_end_template)

        # End XML file.
        output.write(end_template.format(image_data=im_data_template.format(filename=filename, folder=folder)))

if __name__ == '__main__':
    main()