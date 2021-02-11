import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons


def plot_colormap_interactive(x, y, z, xlabel=None, ylabel=None, zlabel=None, additional_cmaps = []):
    '''
    Plots x, y, z, possibly with labels.
    Built-in list of colormaps to choose from: viridis, magma, RdBu_r.
    More cmap options can be added.
    Returns the fig, ax created by plt.subplots() and the image by pcolormesh().
    '''
    
    fig, ax = plt.subplots()
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
    slider_gamma = Slider(ax_gamma, r'$\gamma$', 0, 3, 1)

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
                                        ['viridis', 'magma', 'RdBu_r'] + additional_cmaps,\
                                         activecolor='#206ba9')

    def _colorfunc(label):
        heat_map.set_cmap(buttons_colorscheme.value_selected)
        fig.canvas.draw_idle()

    buttons_colorscheme.on_clicked(_colorfunc)

    return fig, ax, heat_map


# to do: draw and label things better
def plot_linecut_viewers(ax_2d, x, y, z, marker_color='white'):
    """
    Takes the fig, ax, mesh of a heatmap and creates interactive x and y linecut viewers.
    x, y, z need to be supplied and thus are decoupled from actual 2d heatmap data,
    so make sure they agree!
    x, y are vectors and z is a matrix.
    Returns the fig and axes created here.
    """
    fig_1d, axes_1d = plt.subplots(2, 1, figsize=(6,6))
    plt.subplots_adjust(bottom=0.2)

    hline = ax_2d.axhline(y=y[0], color=marker_color)
    vline = ax_2d.axvline(x=x[0], color=marker_color)

    line_xcut, = axes_1d[0].plot(x, z[0, :])
    axes_1d[0].set_xlabel('x')
    line_ycut, = axes_1d[1].plot(y, z[:, 0])
    axes_1d[1].set_xlabel('y')

    ax_yslider = plt.axes([0.2, 0.1, 0.7, 0.03])
    ax_xslider = plt.axes([0.2, 0.05, 0.7, 0.03])
    slider_y = Slider(ax_yslider, 'y index', 0, len(y), valinit=0, valstep=1)
    slider_x = Slider(ax_xslider, 'x index', 0, len(x), valinit=0, valstep=1)

    def _update(val):
        x_ind = slider_x.val
        y_ind = slider_y.val
        data_xcut = z[int(y_ind), :]
        data_ycut = z[:, int(x_ind)]

        line_xcut.set_data(x, data_xcut)
        x_span = data_xcut.max() - data_xcut.min()
        axes_1d[0].set_ylim((data_xcut.min() - x_span*0.05, data_xcut.max() + x_span * 0.05))
        line_ycut.set_data(y, data_ycut)
        y_span = data_ycut.max() - data_ycut.min()
        axes_1d[1].set_ylim((data_ycut.min() - y_span*0.05, data_ycut.max() + y_span * 0.05))

        hline.set_ydata([y[int(y_ind)], y[int(y_ind)]])
        vline.set_xdata([x[int(x_ind)], x[int(x_ind)]])
    #     fig.canvas.draw_idle()

    slider_x.on_changed(_update)
    slider_y.on_changed(_update)

    return fig_1d, axes_1d
