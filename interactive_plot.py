import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons


def plot_colormap_interactive(fig, ax, x, y, z, xlabel=None, ylabel=None, zlabel=None):
    '''
    Takes the fig and ax created from plt.subplots() and plots x, y, z
    '''
    
    plt.subplots_adjust(bottom=0.3, left=0.35)

    heat_map = plt.pcolormesh(x, y, z)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
#     ax.margins(x=0)

    ax_min = plt.axes([0.3, 0.1, 0.5, 0.03])
    ax_max = plt.axes([0.3, 0.15, 0.5, 0.03])
    ax_gamma = plt.axes([0.3, 0.05, 0.5, 0.03])
    ax_colorscheme = plt.axes([0.05, 0.05, 0.15, 0.8])

    cb = plt.colorbar(mappable=heat_map, ax=ax, label=zlabel)

    slider_min = Slider(ax_min, 'Min', np.min(z), np.max(z), valinit=np.min(z))
    slider_max = Slider(ax_max, 'Max', np.min(z), np.max(z), valinit=np.max(z))
    slider_gamma = Slider(ax_gamma, '$\gamma$', 0, 3, 1)

    def _update(val):
        val_min = slider_min.val
        val_max = slider_max.val
        val_gamma = slider_gamma.val
        norm = mpl.colors.PowerNorm(val_gamma, vmin=val_min, vmax=val_max)
        heat_map.set_norm(norm)
        cb.update_normal(heat_map)
        fig.canvas.draw_idle()

    slider_min.on_changed(_update)
    slider_max.on_changed(_update)
    slider_gamma.on_changed(_update)

    buttons_colorscheme = RadioButtons(ax_colorscheme, \
                                        ['viridis', 'plasma', 'magma', \
                                        'cividis', 'seismic', 'gist_stern'],\
                                         activecolor='#206ba9')

    def _colorfunc(label):
        heat_map.set_cmap(buttons_colorscheme.value_selected)
        fig.canvas.draw_idle()

    buttons_colorscheme.on_clicked(_colorfunc)
