# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:27:56 2019
@author: Alberto, Guan
"""

import gdspy
import os, copy
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
import design_tools as tools

# def empty_gds_library():
    # keys = [key for key in gdspy.current_library.cell_dict.keys()]
    # [gdspy.current_library.cell_dict.pop(key) for key in keys]
# empty_gds_library()

def generate_2d_array(title, x_labels, y_labels, unit_function, arg_2d_list, dx=40, dy=40, preview=False):
    """
    Generates a 2D array of patterns with unique identifiers, a title, a column of labels on the left and a row of labels on top
    :param title: string, also filename
    :param x_labels: list of strings, each is the label for a column
    :param y_labels: list of strings, each is the label for a row
    :param unit_function: a function that takes each args tuple as parameters and returns a gdspy shape
    :param dx: float, x distance in um
    :param dy: float, y distance in um
    :return: writes a .gds file and launches the gdspy layout viewer
    """

    if len(x_labels) != len(arg_2d_list[0]) or len(y_labels) != len(arg_2d_list):
        print('Check dimensions of arg_2d_list and x/y_labels!')
        raise ValueError

    cell = gdspy.Cell(title)
    tools.make_title(title, cell, y=2 * dy, text_size=20)

    for y_count, arg_1d_list in enumerate(arg_2d_list):
        for x_count, args in enumerate(arg_1d_list):
            x = x_count * dx
            y = -y_count * dy
            wire_cell = []  # will contain a wire and eventually labels, SmartWalls...
            
            # labels
            if y == 0:
                x_label = x_labels[x_count]
                tools.make_x_label(x_label, wire_cell, y=dy, text_size=5)
            if x == 0:
                y_label = y_labels[y_count]
                tools.make_y_label(y_label, wire_cell, x=-dx/2, y=dy/2, text_size=4)

            wire_shape = unit_function(**args)
            if type(wire_shape) == tuple:
                for single_shape in wire_shape:
                    wire_cell.append(single_shape)
            else:
                wire_cell.append(wire_shape)

            if x_count != 0: # left most row can be labeled by the big X labels
                tools.add_label('X{}Y{}'.format(str(x_count), str(y_count)), wire_cell, x=-dx/2, y=dy/2, text_size=2, layer=11)  # adds a number label to the wire_cell
            
            tools.move_list(wire_cell, x, y)  # translates the entire wire_cell
            cell.add(wire_cell)
            
    gdspy.write_gds(os.path.join(pattern_directory, title + '.gds'), cells=[cell], unit=1.0e-6, precision=1e-11)
    # 1um unit, 1e-11 precision allows total cell area up to ~20mm*20mm
    if preview:
        return gdspy.LayoutViewer(cells=[cell])

def generate_hall_bar_array(title, width_list, length_xx_list, dx, dy, orientation=0, leg_type='default'):
    """
    Generates a gdspy cell populated with single nanowires with the specified parameters
    :param title: str. Also filename
    :param width_list: list of floats. In um
    :param length_xx_list: list of floats. In um
    :param dx: float. Distance in um between neighboring structures
    :param dy: float. Distance in um between neighboring structures
    :param orientation: float in degrees. 0deg is horizontal.
    :param leg_type: {'fishbone', 'default'}
    :return: gdspy viewer
    """

    cell = gdspy.Cell(title)
    tools.make_title(title, cell, y=2 * dy, text_size=20)

    for i, w in enumerate(width_list):
        for j, L in enumerate(length_xx_list):
            wire_list = []  # will contain a wire and eventually labels, SmartWalls...
            x = j * dx
            y = - i * dy

            # labels
            if y == 0:
                x_label = 'L_xx=' + str(L)
                tools.make_x_label(x_label, wire_list, y=dy)
            if x == 0:
                y_label = 'w=' + str(int(w * 1000)) + 'nm\n'
                tools.make_y_label(y_label, wire_list, x=-dx)

            length_S_to_junction = 1.5
            length_total = L + 2 * length_S_to_junction
            length_leg = 1

            if j % 2 == 0:
                reservoir_enable = True
            else:
                reservoir_enable = False

            struct = tools.modified_hall_bar(w, length_total=length_total, length_leg=length_leg, length_xx=L,
                                             leg_offset=0, orientation=orientation, leg_type=leg_type,
                                             reservoir_enable=reservoir_enable)

            wire_list.append(struct)

            tools.add_label('X{}Y{}'.format(str(j), str(i)), wire_list, text_size=2, layer=20)  # adds a number label to the wire_cell
            # tools.add_label('W{}L{}X{}Y{}'.format(str(int(w*1000)), str(L), str(j), str(i)), wire_list, text_size=2)  # adds a number label to the wire_cell
            
            tools.move_list(wire_list, x, y)  # translates the entire wire_cell
            cell.add(wire_list)

    gdspy.write_gds(os.path.join(pattern_directory, title + '.gds'), cells=[cell], unit=1.0e-6, precision=1e-11)
    # 1um unit, 1e-11 precision allows total cell area up to ~20mm*20mm
    return gdspy.LayoutViewer(cells=[cell])

def generate_arg_2d_list(raw_dict, x_key, y_key):
    """
    Takes a raw_dict of the form {..., ..., x_key:x_list, ..., y_key:y_list, ...}
    with the location of x_list and y_list specified by x_index and y_index
    Returns a 2D array of dicts [[{...}, {}, ], [{...}, {}, ], ..., [{...}, ]], 
    each dict being a complete set of arguments as the input of a shape generation function,
    with x_key:x_list replaced by x_key:individual_value_in_list
    """
    x_list, y_list = raw_dict[x_key], raw_dict[y_key]
    assert type(x_list) == list and type(y_list) == list, 'Check if the lists of x and y are correctly labeled by x_key and y_key'

    array_2d = []
    for y in y_list:
        new_row = []
        for x in x_list:
            new_item = copy.deepcopy(raw_dict)
            new_row.append(new_item)
            new_row[-1][x_key] = x
            new_row[-1][y_key] = y
        array_2d.append(new_row)
    return array_2d


'''
======== main ==============
'''

with open('CDI/2D_structures_symmetric.yaml', 'r') as f:
    config_data = load(f, Loader=Loader)
try:
    dx = config_data['dx']
except:
    dx = 50
try:
    dy = config_data['dy']
except:
    dy = 40

pattern_directory = config_data['directory']
for a in config_data['arrays']:
    file_name = a['file_name'].format(**a['parameters']) 
    # allows file names such as 'l1={l1}um' as long as an l1 item exists in the list of parameters
    x_key, y_key = a['x_key'], a['y_key']
    x_labels = [a['x_labels'].format(X=x).strip() for x in a['parameters'][x_key]]
    y_labels = [a['y_labels'].format(Y=y).strip() for y in a['parameters'][y_key]]
    arg_2d_list = generate_arg_2d_list(a['parameters'], x_key, y_key)
    if a['enable'] == True:
        generate_2d_array(file_name, x_labels, y_labels, tools.loop_with_sw, arg_2d_list, dx=dx, dy=dy)
