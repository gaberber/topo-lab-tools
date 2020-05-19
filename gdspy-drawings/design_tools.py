# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 12:18:42 2019
@author: Alberto
This file is a library for gdspy.
It contains functions to build nanowire structures, SmartWalls and labels.
"""
import math
import gdspy as gp


def make_title(name, cell, x=0, y=100, text_size=10, layer=0):
    """
    :param name: str
    :param cell: gp.Cell
    :param x: float
    :param y: float
    :param text_size: float
    :param layer: float
    :return: gp.Cell
    """

    # type checks
    if type(name) != str:
        print("ERROR: first argument must be a string\n")
        return cell
    if type(cell) != gp.Cell:
        print("ERROR: second argument must be a gdspy Cell\n")
        return cell
    title = gp.Text(name, text_size, (x, y), layer=layer)
    cell.add(title)
    return cell


def make_x_label(name, struct_list, x=0, y=50, text_size=5, layer=0):
    # type checks
    if type(name) != str:
        print("ERROR: first argument must be a string\n")
        return struct_list
    label = gp.Text(name, text_size, (x, y), layer=layer)
    struct_list.append(label)
    return struct_list


def make_y_label(name, struct_list, x=-50, y=0, text_size=5, layer=0, align_right=True):
    # type checks
    if type(name) != str:
        print("ERROR: first argument must be a string\n")
        return struct_list
    label = gp.Text(name, text_size, (x, y), layer=layer)
    box = label.get_bounding_box()
    if align_right == True:
        x_shift = box[0][0] - box[1][0]
        label.translate(x_shift, 0)
    struct_list.append(label)
    return struct_list


def add_label(label_text, struct_list, x=0, y=15, text_size=3, layer=0):
    number = gp.Text(label_text, text_size, (x, y), layer=layer)
    struct_list.append(number)
    return struct_list


def move_list(lista, x, y):
    for i, struct in enumerate(lista):
        try:
            lista[i] = struct.translate(x, y)
        except:
            print('GDSPY cannot move object: ', type(struct))
    return lista


def single_wire(w, L, layer=1, orientation=0):
    """
    Returns a gdspy shape with one single wire centered at origin
    """
    angle_rad = orientation * math.pi / 180
    p = [(0, 0), (L * math.cos(angle_rad), L * math.sin(angle_rad))]
    wire = gp.PolyPath(p, w, layer=layer)
    return wire


def parallel_wire(w, L, distance=1, layer=1, orientation=0):
    """
    Returns a gdspy shape with two parallel wires centered at origin
    """
    w1 = single_wire(w, L, layer=layer, orientation=0).translate(0, distance / 2)
    w2 = single_wire(w, L, layer=layer, orientation=0).translate(0, -distance / 2)
    angle_rad = orientation * math.pi / 180
    return gp.fast_boolean(w1, w2, 'or', layer=layer).rotate(angle_rad)


def modified_hall_bar(width, length_total=5, length_leg=1, length_xx=2, leg_offset=0, orientation=0,
                      leg_type='default', reservoir_enable=False,  layer=1):
    """
    :param width: nanowire width
    :param length_total: left to right
    :param length_leg: each leg
    :param length_xx: between legs
    :param leg_offset: allows to move upper leg to the left and lower leg to the right a bit
    :param orientation: degrees. 0 = horizontal
    :param leg_type: 'default' or 'fishbone'
    :param reservoir_enable: bool
    :param layer: shapes on 1 by default
    :return: gdspy object
    """
    leg_type_checklist = {'fishbone', 'default'}
    if leg_type not in leg_type_checklist:
        print('Hall bar type must be in ' + str(leg_type_checklist))
        raise ValueError

    p = [(-length_total / 2, 0), (length_total / 2, 0)]
    if leg_type == 'fishbone':
        p1 = [(0, 0), (-length_leg * math.cos(math.pi / 3), length_leg * math.sin(math.pi / 3))]
    else:
        p1 = [(0, 0), (length_leg * math.cos(math.pi / 3), length_leg * math.sin(math.pi / 3))]
    p2 = [(0, 0), (-length_leg * math.cos(math.pi / 3), -length_leg * math.sin(math.pi / 3))]

    b = gp.PolyPath(p, width, layer=layer)
    leg1 = gp.PolyPath(p1, width, layer=layer).translate(-leg_offset / 2, 0)
    leg2 = gp.PolyPath(p2, width, layer=layer).translate(leg_offset / 2, 0)
    slash1 = gp.fast_boolean(leg1, leg2, "or", layer=layer).translate(-length_xx / 2, 0)
    slash2 = gp.fast_boolean(leg1, leg2, "or", layer=layer).translate(length_xx / 2, 0)
    h = gp.fast_boolean(slash1, slash2, "or", layer=layer)
    hb = gp.fast_boolean(h, b, "or", layer=layer)

    reservoir_width = 2
    reservoir_height = 0.5
    if reservoir_enable:
        hb = add_reservoir_to_shape(hb, reservoir_width, reservoir_height,
                                    position=(-(length_total + reservoir_width) / 2, 0))
        hb = add_reservoir_to_shape(hb, reservoir_width, reservoir_height,
                                    position=((length_total + reservoir_width) / 2, 0))

    angle_rad = orientation * math.pi / 180
    return hb.rotate(angle_rad)


def teleportation_loop(width, l1, l2, la=1, orientation=0, loop_type='mirror-symmetric', reservoir_enable=False, layer=10, sw_enable=False, sw_layer=20, sw_distance=1, twin_enable=True):
    """
    :param width: nanowire width
    :param l1: length of shorter (bottom) side 
    :param l2: length of slanted side
    :param la: length of protruding arms for contacts
    :param orientation: degrees. 0deg = horizontal
    :param loop_type: {'mirror-symmetric', 'center-symmetric', 'asymmetric'}
    :param reservoir_enable: not used for now
    :param layer: shapes on 1 by default
    :return: gdspy object
    """

    allowed_types = {'mirror-symmetric', 'center-symmetric', 'asymmetric'}
    if loop_type not in allowed_types:
        print('loop_type must be in ' + allowed_types)
        raise ValueError
    
    l2_proj_x = l2 * math.cos(math.pi / 3)
    l2_proj_y = l2 * math.sin(math.pi / 3)

    if loop_type == 'asymmetric':
        p1 = [(-l2_proj_x - la * math.cos(math.pi/3), l2_proj_y + la * math.sin(math.pi/3)), (0, 0), (l1, 0), (l1 + l2_proj_x, l2_proj_y)]
        p2 = [(-l2_proj_x, l2_proj_y), (l1 + l2_proj_x + la, l2_proj_y)]
    elif loop_type == 'center-symmetric':
        p1 = [(-la, 0), (l1, 0), (l1 + l2_proj_x, l2_proj_y)]
        p2 = [(0, 0), (l2_proj_x, l2_proj_y), (l2_proj_x + l1 + la, l2_proj_y)]
    else:
        p1 = [(-l2_proj_x, l2_proj_y), (0,0), (l1, 0), (l1 + l2_proj_x, l2_proj_y)]
        p2 = [(-l2_proj_x - la, l2_proj_y), (l1 + l2_proj_x + la, l2_proj_y)]

    path1 = gp.FlexPath(p1, width)
    path2 = gp.FlexPath(p2, width)
    loop = gp.fast_boolean(path1, path2, 'or', layer=layer)
    if twin_enable: # only works on mirro-symmetric types
        loop2 = gp.copy(loop)
        loop2.mirror((l1+l2_proj_x+ 1.5*la, 0), (l1+l2_proj_x+ 1.5*la, 1))
        loop = gp.fast_boolean(loop, loop2, 'or', layer=layer)
        loop = gp.fast_boolean(loop, gp.FlexPath([(l1+l2_proj_x+la, l2_proj_y), (l1+l2_proj_x+3*la, l2_proj_y)], width), 'or', layer=layer)

    if reservoir_enable:
        reservoir_width = 2
        reservoir_height = 0.5
        if loop_type == 'asymmetric':
            left_coordinates = (-l2_proj_x - (la + reservoir_width/2) * math.cos(math.pi/3), l2_proj_y + (la + reservoir_width/2) * math.sin(math.pi/3))
            right_coordinates = (l1 + l2_proj_x + la + reservoir_width/2, l2_proj_y)
            loop = add_reservoir_to_shape(loop, reservoir_width, reservoir_height, position=left_coordinates, rotation=-60, layer=layer)
            loop = add_reservoir_to_shape(loop, reservoir_width, reservoir_height, position=right_coordinates, layer=layer)
        elif loop_type == 'center-symmetric':
            left_coordinates = (-la - reservoir_width/2, 0)
            right_coordinates = (l2_proj_x + l1 + la + reservoir_width/2, l2_proj_y)
            loop = add_reservoir_to_shape(loop, reservoir_width, reservoir_height, position=left_coordinates, layer=layer)
            loop = add_reservoir_to_shape(loop, reservoir_width, reservoir_height, position=right_coordinates, layer=layer)
        else:
            left_coordinates = (-l2_proj_x - la - reservoir_width/2, l2_proj_y)
            right_coordinates = (l1 + l2_proj_x + la + reservoir_width/2, l2_proj_y)
            loop = add_reservoir_to_shape(loop, reservoir_width, reservoir_height, position=left_coordinates, layer=layer)
            loop = add_reservoir_to_shape(loop, reservoir_width, reservoir_height, position=right_coordinates, layer=layer)

    wall_polygon = [(0, -sw_distance), (0, -sw_distance-1), (l1, -sw_distance-1), (l1, -sw_distance), (l1+l2_proj_x+la+2, -sw_distance)]
    wall_polygon += straight_to_zigzag_line([(l1+l2_proj_x+la+2, -sw_distance), (l1+l2_proj_x+la+2, -sw_distance-2), (-l2_proj_x-la-2, -sw_distance-2), (-l2_proj_x-la-2, -sw_distance)], 0.3, 30)

    symmetric_wall = gp.Polygon(wall_polygon, layer=sw_layer)
    if twin_enable:
        wall2 = gp.copy(symmetric_wall)
        wall2.mirror((l1+l2_proj_x+ 1.5*la, 0), (l1+l2_proj_x+ 1.5*la, 1))
        symmetric_wall = gp.fast_boolean(symmetric_wall, wall2, 'or', layer=sw_layer)


    if sw_enable:
        return loop, symmetric_wall # no 30deg option with SWs
    else:
        return loop.rotate(orientation * math.pi / 180)


def loop_with_sw(width, l1, l2, d1, d2, la=1, loop_type='high-angle-60deg', layer=10, sw_layer=20, twin_enable=True):
    """
    The above function got too complicated so started a separate loop generating function for those with smart walls
    :param width: nanowire width
    :param l1: length of shorter (bottom) side 
    :param l2: length of slanted side
    :param d1: distance between top arm and the two side SWs
    :param d2: distance between top arm and the island boundary
    :param la: length of protruding arms for contacts. Check the geometries again if you change this value to anything but the default 1um.
    :param loop_type: {'high-angle', 'grazing-angle'}
    :param layer: shapes on 10 by default
    :param sw_layer: SWs on 20 by default
    :param twin_enable: mirror the structure overall and connect with a 3um segment in the middle
    :return: two gdspy objects: first one containing the wire, the second the SWs
    """

    allowed_types = {'high-angle-symmetric', 'high-angle-60deg', 'high-angle-90deg', 'grazing-angle'}
    if loop_type not in allowed_types:
        print('loop_type must be in ' + allowed_types)
        raise ValueError
    
    # convenience variables
    if loop_type == 'high-angle-90deg':
        l2_proj_x = 0
        l2_proj_y = l2
        island_boundary_horizontal_offset = 0.2 
        # islands need to stay inside arm and some distance away from the joints. This is automatically guaranteed in 60deg configurations
    else:
        l2_proj_x = l2 * math.cos(math.pi / 3)
        l2_proj_y = l2 * math.sin(math.pi / 3)
        island_boundary_horizontal_offset = 0.1

    # drawing a unit loop
    if loop_type == 'high-angle-60deg': # this one is special: only twinned structures for now
        la_proj_x = la * math.cos(math.pi / 3)
        la_proj_y = la * math.sin(math.pi / 3)
        ld = 5 # distance between two units
        p1 = [(l2_proj_x+la_proj_x, l2_proj_y+la_proj_y), (0, 0), (l1, 0), (l1+l2_proj_x, l2_proj_y)]
        p2 = [(l1+l2_proj_x+ld, l2_proj_y), (l1+ld, 0), (2*l1+ld, 0), (2*l1+l2_proj_x+ld+la_proj_x, l2_proj_y+la_proj_y)]
        p3 = [(l2_proj_x, l2_proj_y), (2*l1+l2_proj_x+ld, l2_proj_y)]
        path1 = gp.FlexPath(p1, width)
        path2 = gp.FlexPath(p2, width)
        path3 = gp.FlexPath(p3, width)
        loop = gp.boolean(gp.boolean(path1, path2, 'or', layer=layer), path3, 'or', layer=layer)
    else: # 'high-angle-symmetric', 'high-angle-90deg', 'grazing-angle'
        p1 = [(-l2_proj_x, l2_proj_y), (0,0), (l1, 0), (l1 + l2_proj_x, l2_proj_y)]
        p2 = [(-l2_proj_x - la, l2_proj_y), (l1 + l2_proj_x + la, l2_proj_y)]
        path1 = gp.FlexPath(p1, width)
        path2 = gp.FlexPath(p2, width)
        loop = gp.boolean(path1, path2, 'or', layer=layer)

    # attach smart walls
    if loop_type in {'high-angle-symmetric', 'high-angle-90deg'}:
        left_wall = [(-1, l2_proj_y-d2), (0, l2_proj_y-d2), (-math.cos(math.pi/6), l2_proj_y-d2-0.8)]
        left_wall += straight_to_zigzag_line([(-math.cos(math.pi/6), l2_proj_y-d2-0.8), (-3, l2_proj_y-d2-0.8), (-3, l2_proj_y-d2), (-1, l2_proj_y-d2)], 0.2, 30)
        left_wall = gp.Polygon(left_wall, layer=sw_layer)
        left_wall.translate(island_boundary_horizontal_offset, 0)

        right_wall = [(-1, l2_proj_y-d2), (0, l2_proj_y-d2), (-math.cos(math.pi/6), l2_proj_y-d2-0.8)]
        right_wall += straight_to_zigzag_line([(-math.cos(math.pi/6), l2_proj_y-d2-0.8), (-2, l2_proj_y-d2-0.8), (-2, l2_proj_y-d2), (-1, l2_proj_y-d2)], 0.2, 30)
        right_wall = gp.Polygon(right_wall, layer=sw_layer)
        right_wall.mirror((l1/2, 0), (l1/2, 1))
        right_wall.translate(-island_boundary_horizontal_offset, 0)

        middel_wall = gp.Polygon([(-0.2+island_boundary_horizontal_offset, l2_proj_y-d1), (l1+0.2-island_boundary_horizontal_offset, l2_proj_y-d1), (l1+1, l2_proj_y-d1-1), (-1, l2_proj_y-d1-1)], layer=sw_layer)
        
        symmetric_wall = gp.boolean(left_wall, right_wall, 'or', layer=sw_layer)
        symmetric_wall = gp.boolean(symmetric_wall, middel_wall, 'or', layer=sw_layer)
    elif loop_type == 'high-angle-60deg':
        left_wall = [(2*l2_proj_x-0.5, l2_proj_y-d2), (2*l2_proj_x, l2_proj_y-d2), (2*l2_proj_x-0.5, l2_proj_y-d2-0.5)]
        left_wall += straight_to_zigzag_line([(2*l2_proj_x-0.5, l2_proj_y-d2), (-2, l2_proj_y-d2), (-2, l2_proj_y-d2-1), (2*l2_proj_x-0.5, l2_proj_y-d2-0.5)], 0.2, 30)
        left_wall = gp.Polygon(left_wall, layer=sw_layer)

        right_wall = [(l1+0.5, l2_proj_y-d2), (l1, l2_proj_y-d2), (l1+0.5, l2_proj_y-d2-0.5)]
        right_wall += straight_to_zigzag_line([(l1+0.5, l2_proj_y-d2), (l1-2*l2_proj_x+2.5, l2_proj_y-d2), (l1-2*l2_proj_x+2, l2_proj_y-d2-1), (l1+0.5, l2_proj_y-d2-0.5)], 0.2, 30)
        right_wall = gp.Polygon(right_wall, layer=sw_layer)

        middel_wall = gp.Polygon([(2*l2_proj_x-0.2, l2_proj_y-d1), (l1+0.2, l2_proj_y-d1), (l1+1.5, l2_proj_y-d1-3**0.5), (2*l2_proj_x-1.5, l2_proj_y-d1-3**0.5)], layer=sw_layer)

        symmetric_wall = gp.boolean(left_wall, right_wall, 'or', layer=sw_layer)
        symmetric_wall = gp.boolean(symmetric_wall, middel_wall, 'or', layer=sw_layer)

        twin_wall = gp.copy(symmetric_wall)
        twin_wall.mirror((l2_proj_x+l1+0.5*ld, 0), (l2_proj_x+l1+0.5*ld, 1))
        symmetric_wall = gp.boolean(symmetric_wall, twin_wall, 'or', layer=sw_layer)

    else: # loop-type: 'grazing-angle'
        loop.mirror((0,l2_proj_y/2), (1,l2_proj_y/2))
        left_wall = [(0,-d1), (-1, -d1-math.sqrt(3)), (-3, -d1-math.sqrt(3)), (-4, -d1)]
        left_wall = gp.Polygon(left_wall, layer=sw_layer)
        right_wall = gp.copy(left_wall)
        right_wall.mirror((l1/2, 0), (l1/2, 1))
        symmetric_wall = gp.boolean(left_wall, right_wall, 'or', layer=sw_layer)

    # make twin structures to facilitate growth
    if twin_enable and loop_type != 'high-angle-60deg': # only works on mirror-symmetric types
        loop2 = gp.copy(loop)
        loop2.mirror((l1+l2_proj_x+ 2.5*la, 0), (l1+l2_proj_x+ 2.5*la, 1))
        loop = gp.boolean(loop, loop2, 'or', layer=layer)
        if loop_type in {'high-angle-symmetric', 'high-angle-90deg'}:
            loop = gp.boolean(loop, gp.FlexPath([(l1+l2_proj_x+la, l2_proj_y), (l1+l2_proj_x+4*la, l2_proj_y)], width), 'or', layer=layer)
        else: # 'grazing-angle'
            loop = gp.boolean(loop, gp.FlexPath([(l1+l2_proj_x+la, 0), (l1+l2_proj_x+4*la, 0)], width), 'or', layer=layer)

        wall2 = gp.copy(symmetric_wall)
        wall2.mirror((l1+l2_proj_x+ 2.5*la, 0), (l1+l2_proj_x+ 2.5*la, 1))
        symmetric_wall = gp.boolean(symmetric_wall, wall2, 'or', layer=sw_layer)

    shadow = _add_shadow(symmetric_wall)
    return (loop, symmetric_wall, shadow) if symmetric_wall != None else loop # no 30deg option with SWs

def _add_shadow(sw_shape, shadow_vector=[0, 2.3], sw_layer=20, shadow_layer=14):
    p_list = []
    for individual_polygon in sw_shape.polygons:
        for k, point in enumerate(individual_polygon):
            point1 = point
            if k+1 < len(individual_polygon):
                point2 = individual_polygon[k+1] 
            else:
                point2 = individual_polygon[0]
            point3 = (point1[0] + shadow_vector[0], point1[1] + shadow_vector[1])
            point4 = (point2[0] + shadow_vector[0], point2[1] + shadow_vector[1])
            p_list.append(gp.Polygon([point1, point3, point4, point2]))
    
    final_shadow = p_list[0]
    for parallelogram in p_list[1:]:
        final_shadow = gp.boolean(final_shadow, parallelogram, 'or', layer=shadow_layer)
    final_shadow = gp.boolean(final_shadow, sw_shape, 'not', layer=shadow_layer)
    return final_shadow

def add_reservoir_to_shape(shape, x_width=2, y_height=0.5, position=(0, 0), layer=1, rotation=0):
    """
    Takes a shape and adds a rectangle of specified width and height to some position
    :param shape: gp.PolygonSet
    :param x_width: in um
    :param y_height: in um
    :param position: (x, y) center of the rectangle
    :param layer: shapes are on layer 1 by default
    :param rotation: float in degrees
    :return: None
    """

    reservoir = gp.Rectangle((-x_width/2, -y_height/2), (x_width/2, y_height/2), layer=layer).rotate(rotation * math.pi / 180)

    return gp.fast_boolean(shape, reservoir.translate(position[0], position[1]), 'or', layer=layer)


def zigzag_wire(w, L, l_each, angle=30, orientation=0, layer=1):
    """
    makes a zigzag wire with total length ~= L and each segment having length l_each
    :param w: float, width in um
    :param L: float, total length of the straight line before zigzag-ification
    :param l_each: float, chop the straight line into segments each having length l_each
    :param angle: float in degrees, the angle formed between the new zigzag lines and the original straight line
    :param orientation: float in degrees, that of the original straight line
    :param layer: int, 1 by default for shapes
    :return: gdspy shape containing the zigzag line
    """
    if l_each <= 0 or l_each >= L:
        print('Length of each segment must be 0 < l_each < L')
        raise ValueError

    # make segments
    line_coordinates = [(-L / 2, 0)]  # beginning point
    x = -L/2
    even_count = True
    while x < L/2:
        x += l_each
        even_count = not even_count
        if even_count:
            line_coordinates.append((x, 0))
        else:
            line_coordinates.append((x, -l_each * math.sin(angle * math.pi / 180)))
    
    zigzag_line = gp.FlexPath(line_coordinates, w, layer=layer)

    return zigzag_line.rotate(orientation * math.pi / 180)


def star(w, L, layer=1):
    p1 = [(-L * math.cos(math.pi / 3), L * math.sin(math.pi / 3)), (0, 0),
          (-L * math.cos(math.pi / 3), -L * math.sin(math.pi / 3))]
    p2 = [(0, 0), (L, 0)]
    arm1 = gp.PolyPath(p1, w, layer=layer)
    arm2 = gp.PolyPath(p2, w, layer=layer)
    star = gp.fast_boolean(arm1, arm2, "or", layer=layer)
    return star


def barred_hexagon(w, L, layer=1):
    p1 = [(-L, 0), (-L * math.cos(math.pi / 3), L * math.sin(math.pi / 3)),
          (+L * math.cos(math.pi / 3), L * math.sin(math.pi / 3)), (L, 0),
          (L * math.cos(math.pi / 3), -L * math.sin(math.pi / 3)),
          (-L * math.cos(math.pi / 3), -L * math.sin(math.pi / 3)), (-L, 0)]
    p2 = [(-3 * L, 0), (3 * L, 0)]
    hexagon = gp.PolyPath(p1, w, layer=layer)
    bar = gp.PolyPath(p2, w, layer=layer)
    struct = gp.fast_boolean(hexagon, bar, "or", layer=layer)
    return struct


def real_w(w, b=100):
    return w + b


def SW_distance_limits_45(SW_h, real_nw_w, real_nw_h=0.080, angle=45):
    angle_rad = angle * math.pi / 180
    min_d = 0 + real_nw_w / 2
    max_d = SW_h / math.tan(angle_rad) - real_nw_h / math.tan(angle_rad) - real_nw_w / 2
    just_uncovered = SW_h / math.tan(angle_rad) + real_nw_w / 2
    print('SW distance range (nanowire inside the shadow): [{0:.3f}, {1:.3f}]'.format(min_d, max_d))
    print('SW uncovering distance (nanowire just left out): {0:.3f}'.format(just_uncovered))
    return min_d, max_d, just_uncovered


def SW_distance_limits(SW_h, real_nw_w, Al_angle, real_nw_h=0.080, InSb_angle=45):
    if Al_angle >= InSb_angle:
        print('Warning: Al_angle >= InSb_angle')
    Al_angle_rad = Al_angle * math.pi / 180
    InSb_angle_rad = InSb_angle * math.pi / 180
    min_d = SW_h / math.tan(InSb_angle_rad) + real_nw_w / 2
    max_d = SW_h / math.tan(Al_angle_rad) - real_nw_h / math.tan(Al_angle_rad) - real_nw_w / 2
    just_uncovered = SW_h / math.tan(Al_angle_rad) + real_nw_w / 2
    print('SW distance range (nanowire inside the shadow): [{0:.3f}, {1:.3f}]'.format(min_d, max_d))
    print('SW uncovering distance (nanowire just left out): {0:.3f}'.format(just_uncovered))
    return min_d, max_d, just_uncovered


def add_rhombus_SmartWall(struct_list, angle=40, l=2, delta_x=0, delta_y=-0.6, layer=2):
    arg = ((angle / 2) * math.pi) / 180
    p = [(0, 0), (l * math.sin(arg), -l * math.cos(arg)), (0, -2 * l * math.cos(arg)),
         (-l * math.sin(arg), -l * math.cos(arg))]
    SW = gp.Polygon(p, layer=layer)
    SW.translate(delta_x, delta_y)
    struct_list.append(SW)
    return struct_list


def add_tre_SmartWall(struct_list, angle=40, l=2, delta_x=0, delta_y=-0.6, layer=2):
    arg = ((angle / 2) * math.pi) / 180
    p = [(0, 0), (l * math.sin(arg), -l * math.cos(arg)), (0, -2 * l * math.cos(arg)),
         (-l * math.sin(arg), -l * math.cos(arg))]
    SW = gp.Polygon(p, layer=layer)
    SW.translate(delta_x, delta_y)
    struct_list.append(SW)
    return struct_list


def one_zigzag(start_point, end_point, angle=30):
    angle_rad = angle * math.pi / 180
    if end_point[0] < start_point[0] or (end_point[0] == start_point[0] and end_point[1] > start_point[1]):
        t = start_point
        start_point = end_point
        end_point = t
    length = math.sqrt((end_point[0] - start_point[0]) ** 2 + (end_point[1] - start_point[1]) ** 2)
    d = length / 2 / math.cos(angle_rad)
    try:
        line_angle = math.atan((end_point[1] - start_point[1]) / (end_point[0] - start_point[0]))
    except ZeroDivisionError:
        line_angle = - math.pi / 2
    new_angle = line_angle - angle_rad
    middle_point = (start_point[0] + d * math.cos(new_angle), start_point[1] + d * math.sin(new_angle))
    return middle_point


def straight_to_zigzag_line(p, d, angle):
    lines = len(p) - 1
    for i in range(lines):
        length = math.sqrt((p[lines - i - 1][0] - p[lines - i][0]) ** 2 + (p[lines - i - 1][1] - p[lines - i][1]) ** 2)
        N = int((length + d) // (2 * d))
        p1 = []
        for n in range(N):
            p1.append((p[lines - i - 1][0] + (2 * n + 1) * (p[lines - i][0] - p[lines - i - 1][0]) * d / length,
                       p[lines - i - 1][1] + (2 * n + 1) * (p[lines - i][1] - p[lines - i - 1][1]) * d / length))
        for n in range(N - 1):
            p2 = one_zigzag(p1[N - n - 2], p1[N - n - 1], angle)
            p1.insert(N - n - 1, p2)
        p[lines - i:lines - i] = p1
    return p


def y_symmetry(p):
    p1 = []
    for i in range(len(p)):
        p1.append((-p[i][0], p[i][1]))
    return p1


def add_square_SmartWall(struct_list, l=0.6, delta_x=0, delta_y=-0.6, layer=2, zigzag=False, zigzag_d=0.3,
                         zigzag_angle=30):
    p = [(0, 0), (l, 0), (l, -l), (0, -l), (0, 0)]
    if zigzag == True:
        p = straight_to_zigzag_line(p, zigzag_d, zigzag_angle)
    SW = gp.Polygon(p, layer=layer)
    SW.translate(delta_x, delta_y)
    struct_list.append(SW)
    return struct_list


def add_tree_SmartWall(struct_list, angle=60, l=2, d=0.4, delta_x=0, delta_y=-0.6, layer=2):
    angle_rad = angle * math.pi / 180
    p1 = [(0, 0), (l * math.sin(angle_rad / 2), -l * math.cos(angle_rad / 2))]
    p1 = straight_to_zigzag_line(p1, d, angle / 2)
    p2 = y_symmetry(p1)
    p2 = p2[::-1]
    p = p1 + p2
    SW = gp.Polygon(p, layer=layer)
    SW.translate(delta_x, delta_y)
    struct_list.append(SW)
    return struct_list


def add_double_SmartWall(struct_list, l=5.5, delta_x=0, delta_y=-0.6, relative_y=-0.8, x_superposition=0.15,
                         min_thickness=0.4, layer=2, zigzag=False, zigzag_d=0.3, zigzag_angle=30):
    if -relative_y < min_thickness + 0.3:
        min_thickness = max(0.2, -relative_y - 0.3)
    p1 = [(0, 0), (0, -min_thickness), (-2, -1.5), (-l, -1.5), (-l, 0), (0, 0)]
    p2 = [(-x_superposition, 0), (l, 0), (l, -2), (-0.5, -2), (-0.5, -1.6), (-x_superposition, 0)]
    if zigzag == True:
        p1 = straight_to_zigzag_line(p1, zigzag_d, zigzag_angle)
        p2 = straight_to_zigzag_line(p2, zigzag_d, zigzag_angle)
    SW1 = gp.Polygon(p1, layer=layer)
    SW1.translate(delta_x, delta_y)
    SW2 = gp.Polygon(p2, layer=layer)
    SW2.translate(delta_x, delta_y + relative_y)
    SW = gp.fast_boolean(SW1, SW2, "or", layer=layer)
    struct_list.append(SW)
    return struct_list


def add_triple_SmartWall(struct_list, l=5.5, delta_x=0, delta_y=-0.6, l_island=3, relative_y=-0.8, x_superposition=0.15,
                         min_thickness=0.4, layer=2, zigzag=False, zigzag_d=0.3, zigzag_angle=30):
    if -relative_y < min_thickness + 0.3:
        min_thickness = max(0.2, -relative_y - 0.3)
    p1 = [(-l_island / 2, 0), (-l_island / 2, -min_thickness), (-l_island / 2 - 2, -1.5), (-l, -1.5), (-l, 0),
          (-l_island / 2, 0)]
    p2 = [(-x_superposition - l_island / 2, 0), (x_superposition + l_island / 2, 0),
          (x_superposition + l_island / 2 + 0.5, -1.6), (x_superposition + l_island / 2 + 0.5, -2),
          (-x_superposition - l_island / 2 - 0.5, -2), (-x_superposition - l_island / 2 - 0.5, -1.6),
          (- x_superposition - l_island / 2, 0)]
    if zigzag == True:
        p1 = straight_to_zigzag_line(p1, zigzag_d, zigzag_angle)
        p2 = straight_to_zigzag_line(p2, zigzag_d, zigzag_angle)
    p3 = y_symmetry(p1)
    SW1 = gp.Polygon(p1, layer=layer)
    SW1.translate(delta_x, delta_y)
    SW2 = gp.Polygon(p2, layer=layer)
    SW2.translate(delta_x, delta_y + relative_y)
    SW3 = gp.Polygon(p3, layer=layer)
    SW3.translate(delta_x, delta_y)
    SW = gp.fast_boolean(SW1, SW2, "or", layer=layer)
    SW = gp.fast_boolean(SW, SW3, "or", layer=layer)
    struct_list.append(SW)
    return struct_list


def t_wire(w, L, l=1, up=False, layer=1):
    p1 = [(-L / 2, 0), (L / 2, 0)]
    if up == False:
        p2 = [(0, 0), (-l * math.cos(math.pi / 3), -l * math.sin(math.pi / 3))]
    else:
        p2 = [(0, 0), (l * math.cos(math.pi / 3), l * math.sin(math.pi / 3))]
    arm1 = gp.PolyPath(p1, w, layer=layer)
    arm2 = gp.PolyPath(p2, w, layer=layer)
    t = gp.fast_boolean(arm1, arm2, "or", layer=layer)
    return t


def m_wire(w, L, l=0.5, d=0.5, up=False, layer=1):
    p1 = [(-L / 2, 0), (L / 2, 0)]
    if up == False:
        p2 = [(0, 0), (-l * math.cos(math.pi / 3), -l * math.sin(math.pi / 3))]
    else:
        p2 = [(0, 0), (l * math.cos(math.pi / 3), l * math.sin(math.pi / 3))]
    arm1 = gp.PolyPath(p1, w, layer=layer)
    arm2 = gp.PolyPath(p2, w, layer=layer)
    arm3 = gp.PolyPath(p2, w, layer=layer).translate(-d, 0)
    arm4 = gp.PolyPath(p2, w, layer=layer).translate(d, 0)
    t = gp.fast_boolean(arm1, arm2, "or", layer=layer)
    n = gp.fast_boolean(arm3, arm4, "or", layer=layer)
    m = gp.fast_boolean(t, n, "or", layer=layer)
    return m


def marker_loc(marker, struct=None, loc=2, dx=0, dy=1):
    struct_box = [[0, 0], [0, 0]]
    if struct != None:
        struct_box = struct.get_bounding_box()
    marker_box = marker.get_bounding_box()
    delta_x = 0
    delta_y = 0
    if loc == 1:
        # it's top right
        delta_x = struct_box[1][0] - marker_box[0][0] + dx
        delta_y = struct_box[1][1] - marker_box[0][1] + dy
    if loc == 2:
        # it's top left
        delta_x = struct_box[0][0] - marker_box[1][0] - dx
        delta_y = struct_box[1][1] - marker_box[0][1] + dy
    if loc == 3:
        # it's bottom left
        delta_x = struct_box[0][0] - marker_box[1][0] - dx
        delta_y = struct_box[0][1] - marker_box[1][1] - dy
    if loc == 4:
        # it's bottom right
        delta_x = struct_box[1][0] - marker_box[0][0] + dx
        delta_y = struct_box[0][1] - marker_box[1][1] - dy
    marker.translate(delta_x, delta_y)
    return marker


def armchair_zigzag_wire(w, L, l_each, angle=60, orientation=0, layer=1):
    """
    makes a zigzag wire with total length ~= L and each segment having length l_each
    :param w: float, width in um
    :param L: float, total length of the straight line before zigzag-ification
    :param l_each: float, chop the straight line into segments each having length l_each
    :param angle: float in degrees, the angle formed between the new zigzag lines and the original straight line
    :param orientation: float in degrees, that of the original straight line
    :param layer: int, 1 by default for shapes
    :return: gdspy shape containing the zigzag line
    """
    if l_each <= 0 or l_each >= L:
        print('Length of each segment must be 0 < l_each < L')
        raise ValueError

    # make segments
    line_coordinates = [(-L / 2, 0)]  # beginning point
    x = -L/2
    y = 0
    state = 0
    while x < L/2:
        if state%2 == 0:
            x += l_each
            line_coordinates.append((x, y))
        elif state%4 == 1:
            x += l_each*math.cos(angle * math.pi / 180)
            y += l_each*math.sin(angle * math.pi / 180)
            line_coordinates.append((x, y))
        else:
            x += l_each*math.cos(angle * math.pi / 180)
            y -= l_each*math.sin(angle * math.pi / 180)
            line_coordinates.append((x, y))
        state = state + 1
    
    zigzag_line = gp.FlexPath(line_coordinates, w, layer=layer)

    return zigzag_line.rotate(orientation * math.pi / 180)

def A_marker(number_of_squares=4, d=1, layer=3):
    square = gp.Rectangle((0, 0), (d, -d), layer=layer)
    marker = gp.Rectangle((0, 0), (d, -d), layer=layer)
    for _ in range(1, number_of_squares):
        square.translate(d, -d)
        marker = gp.fast_boolean(marker, square, "or", layer=layer)
    return marker


def B_marker_bits(n, d=1, layer=4):
    """
    Adds a unique bit-marker according to the following powers of two:
    * 0 1 2
      *   3
        * 4
          *
    If more powers are needed:
    * 0 1 2 8
      *   3 7
        * 4 6
          * 5
            *
    Where the asterisks stays for the A_marker
    """
    if n <= 0:
        print('This funciton works with positive integers only.')
        raise ValueError
    if n > 512:
        print('Hey man, how many markers are you planning to put? This code works only up to 512 markers per row or column.')
        raise ValueError

    square = gp.Rectangle((0, 0), (d, -d), layer=layer).translate(d, 0)
    q = n
    # 2^0
    marker=None
    if q % 2 == 1:
        marker = gp.fast_boolean(square, marker, "or", layer=layer)
    # 2^1
    square.translate(d, 0)
    q = q//2
    if q % 2 == 1:
        marker = gp.fast_boolean(square, marker, "or", layer=layer)
    # 2^2
    square.translate(d, 0)
    q = q//2
    if q % 2 == 1:
        marker = gp.fast_boolean(square, marker, "or", layer=layer)
    # 2^3
    square.translate(0, -d)
    q = q//2
    if q % 2 == 1:
        marker = gp.fast_boolean(square, marker, "or", layer=layer)
    # 2^4
    square.translate(0, -d)
    q = q//2
    if q % 2 == 1:
        marker = gp.fast_boolean(square, marker, "or", layer=layer)
    # 2^5
    square.translate(d, -d)
    q = q//2
    if q % 2 == 1:
        marker = gp.fast_boolean(square, marker, "or", layer=layer)
    # 2^6
    square.translate(0, d)
    q = q//2
    if q % 2 == 1:
        marker = gp.fast_boolean(square, marker, "or", layer=layer)
    # 2^7
    square.translate(0, d)
    q = q//2
    if q % 2 == 1:
        marker = gp.fast_boolean(square, marker, "or", layer=layer)
    # 2^8
    square.translate(0, d)
    q = q//2
    if q % 2 == 1:
        marker = gp.fast_boolean(square, marker, "or", layer=layer)

    return marker


def add_marker(struct_list, d=1, struct=None, loc=2, dx=2, dy=2, number_of_squares=3, layer=3):
    marker = A_marker(number_of_squares=number_of_squares, d=d, layer=layer)
    marker = marker_loc(marker, struct=struct, loc=loc, dx=dx, dy=dy)
    struct_list.append(marker)
    return struct_list


def bitmarker(nx, ny, d=1, layer=3):
    nA = 4
    if (nx >= 32) or (ny >= 32):
        nA = 5
    marker = A_marker(number_of_squares=nA, d=d, layer=layer)
    if nx > 0:
        markerBx = B_marker_bits(nx, d=d, layer=layer)
        marker = gp.fast_boolean(marker, markerBx, "or", layer=layer)
    if ny > 0:
        markerBy = B_marker_bits(ny, d=d, layer=layer).mirror((d,-d))
        marker = gp.fast_boolean(marker, markerBy, "or", layer=layer)
    return marker


def add_bitmarker(struct_list, nx, ny, d=1, struct=None, loc=2, dx=2, dy=2, number_of_squares=3, layer=3):
    marker = bitmarker(nx, ny, d=d, layer=layer)
    marker = marker_loc(marker, struct=struct, loc=loc, dx=dx, dy=dy)
    struct_list.append(marker)
    return struct_list