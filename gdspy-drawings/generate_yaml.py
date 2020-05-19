# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 2020
@author: Di
"""
# inputs
Directory = 'E:/Documents/patterns'
fn = '2D_structures_symmetric.yaml' # yaml file name

# define gds file names
# File_names = ['loop-90deg-l1={l1}um-width={width}um', \
    # 'loop-60deg-l1={l1}um-width={width}um', \
    #     'loop-symmetric-l1={l1}um-width={width}um']
File_names = ['loop-symmetric-l1={l1}um-width={width}um']  # seperate each type in different yaml files - to avoid overflow
Enable_file_names = [1, 1, 1] # 1: enable: True; 0: enable: False. Same size as File_names
Twin_file_names = [1, 1, 1]

# define corresponding types
Loop_types = []
for i, File_name in enumerate(File_names):
    if File_name == 'loop-90deg-l1={l1}um-width={width}um':
        Loop_types.extend(['high-angle-90deg'])
    elif File_name == 'loop-60deg-l1={l1}um-width={width}um':
        Loop_types.extend(['high-angle-60deg'])
    elif File_name == 'loop-symmetric-l1={l1}um-width={width}um':
        Loop_types.extend(['high-angle-symmetric'])

# parameters staying constant in a field
width_settings = [0.02, 0.03, 0.04, 0.05]
Enable_width = [1, 1, 1, 1] # 1: enable: True; 0: enable: False. Same size as width_settings
Twin_width = [1, 1, 1, 1]

l1_settings = [1, 1.3, 1.6, 2]
Enable_l1 = [1, 1, 1, 1] # 1: enable: True; 0: enable: False. Same size as l1_settings
Twin_l1 = [1, 1, 1, 1]

# parameters varying in a field
l2 = '[0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.8]'
d1 = '[2.1, 2.2, 2.3, 2.35, 2.4, 2.5, 2.6]'

# parameters always staying constant
d2 = 1.2
la = 1

layer = 10
sw_layer = 20

'''
======== main ==============
'''

ind = '    '

file = open(fn, 'w+')

file.write('directory: ' + Directory + '\n')
file.write('arrays:\n')

for i1, File_name in enumerate(File_names):
    for i2, width in enumerate(width_settings):
        for i3, l1 in enumerate(l1_settings):
            file.write('-\n')
            file.write(ind + 'file_name: ' + File_names[i1] + '\n')
            if Enable_file_names[i1] == 1 and Enable_width[i2] == 1 and Enable_l1[i3] == 1:
                file.write(ind + 'enable: True\n')
            else:
                file.write(ind + 'enable: False\n')
            file.write(ind + 'y_labels: |\n' + 2*ind + 'd1\n' + 2*ind + '{Y}\n')
            file.write(ind + 'x_labels: l2={X}um\n')
            file.write(ind + 'x_key: l2\n')
            file.write(ind + 'y_key: d1\n')
            file.write(ind + 'parameters:\n')
            file.write(ind + ind + 'width: ' + str(width) + '\n')
            file.write(ind + ind + 'l1: ' + str(l1) + '\n')
            file.write(ind + ind + 'l2: ' + l2 + '\n')
            file.write(ind + ind + 'd1: ' + d1 + '\n')
            file.write(ind + ind + 'd2: ' + str(d2) + '\n')
            file.write(ind + ind + 'la: ' + str(la) + '\n')
            file.write(ind + ind + 'loop_type: ' + Loop_types[i1] + '\n')
            file.write(ind + ind + 'layer: ' + str(layer) + '\n')
            file.write(ind + ind + 'sw_layer: ' + str(sw_layer) + '\n')
            if Twin_file_names[i1] == 1 and Twin_width[i2] == 1 and Twin_l1[i3] == 1:
                file.write(ind + ind + 'twin_enable: True\n')
            else:
                file.write(ind + ind + 'twin_enable: False\n')
