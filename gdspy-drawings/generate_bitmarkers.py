# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:27:56 2019
@author: Alberto, Guan
"""

import os, string
import gdspy as gp
import design_tools as tools

def bitmarker_array():
    dx = 20
    dy = 20 
    # text_size = 100
    nrows = 32
    ncols = 32

    output_dir = r'/Users/wgz/Downloads'
    title = 'BitmarkerArray'

    cell = gp.Cell(title)

    for x in range(ncols):
        for y in range(nrows):
            # x_label = string.ascii_uppercase[x]
            # y_label = str(y)
            # shape = gp.Text(x_label+y_label, text_size, position=(x*dx, -y*dy), layer=10)
            shape = tools.bitmarker(x, y)
            shape.translate(x*dx, -y*dy)
            cell.add(shape)

    gp.LayoutViewer(cells=[cell])

    gp.write_gds(os.path.join(output_dir, title + '.gds'), cells=[cell], unit=1.0e-6, precision=1e-10)

def RP20_set():
    pass

def bitmarker_2inch():
    dx = 20
    dy = 20 
    # text_size = 100
    nrows = 32
    ncols = 32
    output_dir = r'/Users/wgz/Downloads'
    title = 'BitmarkerArray'
    cell = gp.Cell(title)
    for x in range(ncols):
        for y in range(nrows):
            # x_label = string.ascii_uppercase[x]
            # y_label = str(y)
            # shape = gp.Text(x_label+y_label, text_size, position=(x*dx, -y*dy), layer=10)
            shape = tools.bitmarker(x, y)
            shape.translate(x*dx, -y*dy)
            cell.add(shape)
    gp.LayoutViewer(cells=[cell])
    gp.write_gds(os.path.join(output_dir, title + '.gds'), cells=[cell], unit=1.0e-6, 
    precision=1e-10)


