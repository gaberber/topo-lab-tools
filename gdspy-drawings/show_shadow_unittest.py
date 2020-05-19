# -*- coding: utf-8 -*-
"""
@author: Guan
"""

import gdspy as gp
import design_tools as tools


tested_loop, tested_wall, tested_shadow = tools.loop_with_sw(width=0.04, l1=1.3, l2=0.8, d1=2.3, d2=1.2, la=1, loop_type='high-angle-60deg', layer=10, sw_layer=20, twin_enable=True)

c = gp.Cell('shadow_test')
print(type(tested_wall))
c.add(tested_loop)
c.add(tested_wall)
c.add(tested_shadow)

gp.LayoutViewer(cells=[c])
