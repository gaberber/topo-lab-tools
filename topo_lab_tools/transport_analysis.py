import os, re
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import pandas as pd
import copy
from qcodes import load_by_run_spec, initialise_or_create_database_at, load_by_id
import qcodes as qc

mpl.rcParams["text.usetex"] = False
#mpl.rcParams["font.family"] = "Helvetica"
mpl.rcParams["text.latex.preamble"] = "\\usepackage{amsmath}"
mpl.rcParams["figure.max_open_warning"]
mpl.rcParams["axes.labelsize"] = 16
mpl.rcParams["xtick.labelsize"] = 16
mpl.rcParams["ytick.labelsize"] = 16
mpl.rcParams["font.size"] = 15
mpl.rcParams["figure.figsize"] = [12,5]
mpl.rcParams["ytick.major.size"] = 8
mpl.rcParams["ytick.major.width"] = 1.5
mpl.rcParams["ytick.minor.size"] = 4
mpl.rcParams["ytick.minor.width"] = 1
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams["xtick.major.size"] = 8
mpl.rcParams["xtick.major.width"] = 1.5
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["xtick.minor.size"] = 4
mpl.rcParams["xtick.minor.width"] = 1
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["savefig.pad_inches"] = 0.05
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["legend.fontsize"] = 14


"""
---------------------------------------------------------------------------------------------------
                                Generic Functions and Data Objects
---------------------------------------------------------------------------------------------------
"""
def guess_2D_dims(inner_axis, outer_axis, data_2d, rescale_xy = True):
    '''
    Takes X, Y and Z from load_by_id() as inputs.
    Returns reshaped X, Y, Z and the dimensions guessed.
    With incomplete sweeps the unifinished row is discarded.
    '''

    dims_guess = [0, 0]
    dims_guess[1] = sum([outer_axis[0] == kk for kk in outer_axis])
    dims_guess[0] = int(np.floor(len(outer_axis) / dims_guess[1]))
    total_num = dims_guess[1] * dims_guess[0] # for incomplete sweeps where the unifinished row is discarded

    if rescale_xy:
        return inner_axis[:total_num].reshape(dims_guess)[0,:], \
                         outer_axis[:total_num].reshape(dims_guess)[:,0], \
                         data_2d[:total_num].reshape(dims_guess), \
                         dims_guess
    else:
        return data_2d[:total_num].reshape(dims_guess), dims_guess


def batch_guess_2D_dims(inner_axis, outer_axis, data_2d_list):
    '''
    Runs guess_2D_dims on each element of the data_2d_list
    Returns reshaped reshaped inner, outer axes and a list of matrices
    '''
    new_dlist = []
    for d in data_2d_list:
        new_d, dims = guess_2D_dims(inner_axis, outer_axis, d, rescale_xy=False)
        new_dlist.append(new_d)
    
    new_inner_axis, new_outer_axis, _, dims = guess_2D_dims(inner_axis, outer_axis, d, rescale_xy=True)

    return new_inner_axis, new_outer_axis, new_dlist, dims


class Dataset_3T():
    """
    A data-storage-format-agnostic class storing and performing most common operations
    on full conductance matrix measurements.
    """
    def __init__(self, VL, VR, Is_VL, Is_VR, gs_VL, gs_VR, y_param, y_param_label):
        """    
        Inputs are np.array's in proper shapes or lists of them.
        VL, VR: DC bias on left and right. Unit: V
        Is_VL, Is_VR: each contains IL and IR wrt respective bias sweeps. Unit: A
        gs_VL, gs_VR: each contains in order: gLL, gLR, gRL, gRR wrt respective biases. Unit: G0
        y_param, y_param_label: meaning-agnostic y axis.
        Unpacks all inputs into individual np.array's.
        """
        self.VL, self.VR = VL, VR
        self.IL_VL, self.IR_VL = Is_VL
        self.IL_VR, self.IR_VR = Is_VR
        self.gLL_VL, self.gLR_VL, self.gRL_VL, self.gRR_VL = gs_VL
        self.gLL_VR, self.gLR_VR, self.gRL_VR, self.gRR_VR = gs_VR
        self.y_param, self.y_param_label = y_param, y_param_label
        
    def plot_cond_matrix(self, corr = False, sym_cmap=True):
        """
        Creates a figure with four panels plotting the conductance matrix.
        Stores the fig and axes as attributes of this instance for external use.
        Enabling corr option will make it plot gLL_VL_corr etc.
        """
        def _make_cbar(ax, mesh, label):
            cbar = plt.colorbar(mesh, ax=ax)
            cbar.set_label(label, rotation=90)
            return cbar
        
        self.fig, self.axes = plt.subplots(2, 2, figsize = (8, 8))
        self.heatmaps = {}
        
        if corr:
            gLL, gLR, gRL, gRR = self.gLL_VL_corr, self.gLR_VR_corr, self.gRL_VL_corr, self.gRR_VR_corr
            VL, VR = self.VL_corr, self.VR_corr
        else:
            gLL, gLR, gRL, gRR = self.gLL_VL, self.gLR_VR, self.gRL_VL, self.gRR_VR
            VL, VR = self.VL, self.VR

        ax = self.axes[0,0]
        mesh = ax.pcolormesh(VL*1e3, self.y_param, gLL, cmap='magma')
        ax.set_xlabel('$V_L$ (mV)')
        ax.set_title('$g_{LL}$')
        cbar = _make_cbar(ax, mesh, '$g_{LL}$ ($G_0$)')
        self.heatmaps['LL'] = mesh, cbar

        ax = self.axes[0,1]
        mesh = ax.pcolormesh(VR*1e3, self.y_param, gLR*1e3, cmap='RdBu_r')
        ax.set_xlabel('$V_R$ (mV)')
        ax.set_title('$g_{LR}$')
        cbar = _make_cbar(ax, mesh, '$g_{LR}$ ($10^{-3} G_0$)')
        self.heatmaps['LR'] = mesh, cbar
        if sym_cmap:
            vmax = np.abs(gLR*1e3).max()
            vmin = -vmax
            self.update_minmax('LR', vmin, vmax)

        ax = self.axes[1,0]
        mesh = ax.pcolormesh(VL*1e3, self.y_param, gRL*1e3, cmap='RdBu_r')
        ax.set_xlabel('$V_L$ (mV)')
        ax.set_title('$g_{RL}$')
        cbar = _make_cbar(ax, mesh, '$g_{RL}$ ($10^{-3} G_0$)')
        self.heatmaps['RL'] = mesh, cbar
        if sym_cmap:
            vmax = np.abs(gRL*1e3).max()
            vmin = -vmax
            self.update_minmax('RL', vmin, vmax)

        ax = self.axes[1,1]
        mesh = ax.pcolormesh(VR*1e3, self.y_param, gRR, cmap='magma')
        ax.set_xlabel('$V_R$ (mV)')
        ax.set_title('$g_{RR}$')
        cbar = _make_cbar(ax, mesh, '$g_{RR}$ ($G_0$)')
        self.heatmaps['RR'] = mesh, cbar

        for ax in self.axes.flatten():
            ax.set_ylabel(self.y_param_label)
        self.fig.tight_layout()
        
    def update_minmax(self, label, vmin, vmax, gamma=None):
        """
        label: in {'LL', 'LR', 'RL', 'RR'}
        Updates the corresponding panel directly
        """
        if gamma == None:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            assert gamma > 0, 'Check input: gamma needs to be a positive number'
            norm = mpl.colors.PowerNorm(gamma, vmin=vmin, vmax=vmax)
        mesh, cbar = self.heatmaps[label]
        mesh.set_norm(norm)
        cbar.update_normal(mesh)
        
    def correct_dcac(self, RL, RR, RS):
        """
        RL, RR, RS are the total resistances in Ohm on the left, right and middle side, respectively.
        Creates self.VR_corr, self.VL_corr, self.gLL_VL_corr etc etc.
        """
        # correcting DC bias 
        self.VR_corr = (self.VR - self.IR_VR * RR - (self.IL_VR + self.IR_VR) * RS)
        self.VL_corr = (self.VL - self.IL_VL * RL - (self.IL_VL + self.IR_VL) * RS)

        # correcting AC signal
        G0=7.748091191125409e-05
        R_mat = np.array([[RL+RS, RS], [RS, RR+RS]])
        self.gLL_VL_corr, self.gLR_VL_corr, self.gRL_VL_corr, self.gRR_VL_corr = [np.empty_like(self.gLL_VL) for k in range(4)]
        self.gLL_VR_corr, self.gLR_VR_corr, self.gRL_VR_corr, self.gRR_VR_corr = [np.empty_like(self.gLL_VL) for k in range(4)]
        for k_row in range(self.gLL_VL.shape[0]):
            for k_col in range(self.gLL_VL.shape[1]):
                Gmat_VL = G0 * np.array([[self.gLL_VL[k_row, k_col], self.gLR_VL[k_row, k_col]], \
                                  [self.gRL_VL[k_row, k_col], self.gRR_VL[k_row, k_col]]])
                V_jacobian_VL = np.eye(2) - R_mat.dot(Gmat_VL)
                Gmat_VL_corr = Gmat_VL.dot(np.linalg.inv(V_jacobian_VL))
                self.gLL_VL_corr[k_row, k_col], self.gLR_VL_corr[k_row, k_col], self.gRL_VL_corr[k_row, k_col], self.gRR_VL_corr[k_row, k_col] = Gmat_VL_corr.flatten()/G0
                
                Gmat_VR = G0 * np.array([[self.gLL_VR[k_row, k_col], self.gLR_VR[k_row, k_col]], \
                                  [self.gRL_VR[k_row, k_col], self.gRR_VR[k_row, k_col]]])
                V_jacobian_VR = np.eye(2) - R_mat.dot(Gmat_VR)
                Gmat_VR_corr = Gmat_VR.dot(np.linalg.inv(V_jacobian_VR))
                self.gLL_VR_corr[k_row, k_col], self.gLR_VR_corr[k_row, k_col], self.gRL_VR_corr[k_row, k_col], self.gRR_VR_corr[k_row, k_col] = Gmat_VR_corr.flatten()/G0


"""
---------------------------------------------------------------------------------------------------
                                    Spyview Data Objects
---------------------------------------------------------------------------------------------------
"""
def import_spyview_dat(data_dir, filename):
    """
    Returns a np.array in the same shape as the raw .dat file
    """
    with open(os.path.join(data_dir, filename)) as f:
        dat = np.loadtxt(f)
    return dat


class Dataset_2d():
    def __init__(self, x_init, y_init, z_init_list):
        """
        X and Y are arrays. Z is a a list of arrays. 
        All arrays should have the same length.
        X and Y will be then stored as a single 1D array while each element of Z as a matrix.
        """
        x, y, z_list, dims = batch_guess_2D_dims(x_init, y_init, z_init_list)
        self.x = x
        self.y = y
        self.z_list = z_list
        self.dims = dims
        self.create_labels('', '', [''] * len(z_list)) # fill with empty labels first
    
    def create_labels(self, xlabel, ylabel, zlabel_list):
        """
        Creates xlabel, ylabel, zlabel_list
        """
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel_list = zlabel_list
        
    def add_z(self, z, zlabel=''):
        """
        Adds a matrix to z_list along with its label
        """
        self.z_list.append(z)
        self.zlabel_list.append(zlabel)
        
    def plot_all_z(self):
        """
        You probably want to `%matplotlib inline` first
        """
        for k, z in enumerate(self.z_list):
            fig, ax = plt.subplots(1, 1)
            mesh = ax.pcolormesh(self.x, self.y, z)
            fig.colorbar(mesh, ax=ax, label=self.zlabel_list[k])
            ax.set_xlabel(self.xlabel)
            ax.set_ylabel(self.ylabel)


class Dataset_2d_spyview(Dataset_2d):
    def __init__(self, data_dir, filename):
        """
        Reads both .dat and .meta.txt
        Creates a dataset with ALL columns stored in z_list using the metadata in txt.
        No unit conversions — this class is unaware of the meaning of the data.
        To do that, inherit this class and extract values from the general x, y, z_list
        while expressing them in right units.
        """
        self.filename = filename
        dat = import_spyview_dat(data_dir, filename)
        
        x_init = dat[:, 0]
        y_init = dat[:, 1]
        z_init_list = [dat[:, k] for k in range(3, dat.shape[1])]
        Dataset_2d.__init__(self, x_init, y_init, z_init_list)
        
        with open(os.path.join(data_dir, filename[:-3] + 'meta.txt')) as f:
            metadata = f.readlines()
        col_labels = [metadata[k].strip() for k in range(13, len(metadata), 2)]
        
        Dataset_2d.create_labels(self, metadata[3].strip(), metadata[7].strip(), col_labels)

    def plot_columns(self, columns_list):
        """
        columns_list: a list of columns (as shown in Spyview, so starting from 4) to plot
        Returns the fig and axs created here.
        """
        fig, axs = plt.subplots(len(columns_list), 1, figsize=(6, 2*len(columns_list)))
        for ax, col in zip(axs, columns_list):
            if self.dims[0] == 1:
                ax.plot(self.x, self.z_list[col - 4][0,:])
                ax.set_ylabel(self.zlabel_list[col - 4])
                ax.set_xlabel(self.xlabel)
            else:
                mesh = ax.pcolormesh(self.x, self.y, self.z_list[col-4], cmap='viridis')
                ax.set_ylabel(self.ylabel)
                cbar = plt.colorbar(mesh, ax = ax, label=self.zlabel_list[col - 4])
        axs[0].set_title(self.filename)
        fig.tight_layout()
        return fig, axs


class Dataset_BiasVsGate_spyview(Dataset_2d_spyview):
    def __init__(self, data_dir, filename, current_col = 4, conductance_col = 6):
        """
        X: bias. Y: gate.
        Current measured by K1 and conductance by L1 (columns 4 and 6 in spyview by default).
        Parses the metadata txt to interpret the gate and bias amplifications.
        All units in SI.
        """
        Dataset_2d_spyview.__init__(self, data_dir, filename)
        try:
            self.V_range = float(re.findall('[0-9]+mV/V', self.xlabel)[0][:-4]) * 1e-3
            self.gate_amplification = float(re.findall('[0-9]+V/V', self.ylabel)[0][:-3]) * 1e-3
        except:
            print('Check metadata: bias [XX]mV/V and gate [YY]V/V?')
        
        self.Vbias = self.x * 1e-3 * self.V_range # mV -> V then Vsource amp factor
        self.Vg = self.y * self.gate_amplification # assuming gates are always xx V/V
        self.current = self.z_list[current_col - 4]
        self.g = self.z_list[conductance_col - 4]
        
    def correct_bias(self, Rline):
        """
        Needs Rline in Ohm. Creates self.Vbias_cor
        """
        self.Vbias_cor = np.empty_like(self.g)
        for k in range(len(self.Vg)):
            self.Vbias_cor[k, :] = self.Vbias - self.current[k, :] * Rline
        
    def quick_plot_g(self, x_unit = 'mV', bias_cor = True):
        """
        x_unit must be in {'V', 'mV', 'uV'}
        Returns the plotting objects so that you can e.g. mesh.set_norm() later
        """
        unit_mapping = {'V': 1, 'mV': 1e3, 'uV': 1e6}
        x_factor = unit_mapping[x_unit]
        
        fig, ax = plt.subplots()
        if bias_cor:
            mesh = ax.pcolormesh(self.Vbias_cor * x_factor, self.Vg, self.g)
        else:
            mesh = ax.pcolormesh(self.Vbias * x_factor, self.Vg, self.g)
        ax.set_xlabel(f'Bias ({x_unit})')
        ax.set_ylabel('Gate (V)')
        ax.set_title(self.filename)
        fig.colorbar(mesh, ax=ax, label= 'd$I$/d$V$ ($2e^2/h$)')
        
        return fig, ax, mesh


"""
---------------------------------------------------------------------------------------------------
                                    Qcodes Data Objects
---------------------------------------------------------------------------------------------------
"""
        
def init_qcodes(data_dir, db_name):
    qc.config["core"]["db_location"] = data_dir
    assert os.path.isdir(data_dir), "'{}' Is not a Directory".format(data_dir)
    ls= os.listdir(data_dir)
    A = [re.findall(db_name, st.replace(".db", "")) for st in ls]
    i = 0
    for a in A:
        if len(a)>0:
            if a[0]+".db" in ls: 
                i+=1
    assert i>0, "'{0}' not found in '{1}' \n These exist: {2}".format(db_name, data_dir, ls)
    initialise_or_create_database_at("{0}/{1}".format(data_dir, db_name) + ".db" )

class Dataset_qcodes():
    def __init__(self, run_id):
        """
        Reads a qcodes database. automatically renames instruments (keithleys and lockins)
        Arguments:
            run_id: the id associated with the measurement
        Instances: 
            id
            exp_name
            sample_name
            data
            data_labels
            instruments
        Methods:
            subt_line_resist: subtracts line resistance from Vbias, Vac
            diff_cond: computes differential conductance and stores it in self.G
            colormesh: 2D plot of data
            linecut: 1D slice of data
        """
        
        self.dataset = load_by_run_spec(captured_run_id=run_id)
        self.id = run_id
        self.exp_name = self.dataset.exp_name
        self.sample_name = self.dataset.sample_name
        self.data = self.dataset.get_data_as_pandas_dataframe()
        self.data_labels = self.data.keys()
        self.instruments = self.dataset.snapshot['station']["instruments"].keys()
        self._rename_instruments_init()
        self.line_resist_sub = False
        
    def diff_cond(self, column, ac_scaling = 1):
        """computes differential conductance and stores it in self.G
        Works with and without running self.subt_line_resist
        Arguments:
            column: string like 'lockin 1' or 'lockin 2' containing the dI data
            ac_scaling: conversion of amp on lockin to actual unit in Volts
        Creates:
            self.G: the differential conductance in G0 """
        
        assert column in self.data_labels, "column doesn't exist, choose from: {}".format(self.data_labels)
        dI = self.data[column]
        dV = ac_scaling * self.lockins_ac[column]
        if not self.line_resist_sub:
            print("Beware, line resistance is not subtracted yet")
            y,x = self.data[column].index.levels
            self.y = y.values
            self.y_lab = y.name
            self.x = x.values
            self.x_lab = x.name
            self.G = 12906*(dI.values.reshape(len(dI))/dV)
            self.G = self.G.reshape(len(self.y), len(self.x))
            if len(self.y) == len(self.x):
                self.G = self.G[::1, ::-1].T
            
        elif self.line_resist_sub:
            if type(dI.index) == pd.MultiIndex:
                _dI = dI.values.reshape(len(dI))
                L = np.prod(dV.shape)
                if len(_dI) < L:
                    _dI = np.pad(_dI, (0, L-len(_dI)))
                self.G = 12906*_dI.reshape(len(_dI))/dV.reshape(len(_dI))
                try:
                    self.G = self.G.reshape(self.z.shape)
                except (AttributeError,ValueError):
                    x,y = self.data[column].index.levels
                    self.G = self.G.reshape(len(x), len(y))
                    if len(y) == len(x):
                        self.G = self.G[::1, ::-1].T
            else:
                self.G = 12906*dI.values.reshape(len(dI))/dV.reshape(len(dI))
    
    def colormesh(self, column="lockin 1", **pcolorkwargs):
        """2d plot of a measured value
        Arguments:
            column: string like 'lockin 1' or 'keithley 2' 
            pcolorkwargs: [OPT] any argument you would put in plt.pcolor()"""
        if not 'ax' in pcolorkwargs:
            fig, ax = plt.subplots(figsize=(15,6))
        else:
            ax=pcolorkwargs['ax']
            del pcolorkwargs['ax']
        title = self.exp_name + " ({})".format(str(self.id))
        ax.set_title(title)
        if "lockin" in column:
            assert self.G.any(), "Differential Conductance not computed, run diff_cond"
            y = np.repeat(self.y, self.G.shape[1]).reshape(self.G.shape)
            if self.line_resist_sub:
                p = ax.pcolor(self.V_new, y, self.G, cmap = 'inferno', **pcolorkwargs)
            else:
                p = ax.pcolor(self.x, self.y, self.G, cmap = 'inferno', **pcolorkwargs)
            ax.set_xlabel(self.x_lab)
            ax.set_ylabel(self.y_lab)
            cbar = plt.colorbar(p, ax=ax)
            cbar.set_label(r"$dI/dV (G_0)$")
        elif "keithley" in column:
            assert column in self.data_labels, "column doesn't exist, choose from: {}".format(self.data_labels)
            y,x = self.data[column].index.levels
            I_measure = self.data[column].values
            I_measure = I_measure.reshape(len(I_measure))
            L = len(x) * len(y)
            if len(I_measure) < L:
                I_measure = np.pad(I_measure, (0, L-len(I_measure)))
            try:    
                I_measure = I_measure.reshape((len(y), len(x)))
                y = np.repeat(self.y, len(x)).reshape((len(y), len(x)))
            except ValueError:
                I_measure = I_measure.reshape((len(x), len(y)))
                y = np.repeat(self.y, len(y)).reshape((len(x), len(y)))
            if self.line_resist_sub:
                p = ax.pcolor(self.V_new, y, I_measure, cmap = 'inferno')
            else:
                p = ax.pcolor(self.x, y, I_measure, cmap = 'inferno')
            ax.set_xlabel(self.x_lab)
            ax.set_ylabel(self.y_lab)
            cbar = plt.colorbar(p, ax=ax)
            cbar.set_label(column)
            
            
    def linecut(self, column="lockin 1", line=0, digs = 2, **plotkwargs):
        """1d linecut plot of a 2dmap
        Arguments:
            column: string like 'lockin 1' or 'keithley 2' 
            line: [OPT] index of the y-value for which you want a linecut
            digs: [OPT] amount of digits of y-value in legend
            pcolorkwargs: [OPT] any argument you would put in plt.pcolor()"""
        
        title = self.exp_name + " ({})".format(str(self.id))
        plt.title(title)
        if "lockin" in column:
            assert self.G.any(), "Differential Conductance not computed, run diff_cond"
            #y = np.repeat(self.y, self.G.shape[1]).reshape(self.G.shape)
            if self.line_resist_sub:
                p = plt.plot(self.V_new[line,:], self.G[line,:], label = self.y_lab + "= {}".format(round(self.y[line], digs)), **plotkwargs)
            else:
                p = plt.plot(self.x[line, :], self.G[line,:], label = self.y_lab + "= {}".format(round(self.y[line], digs)), **plotkwargs)
            plt.xlabel(self.x_lab)
            plt.ylabel(r"$dI/dV (G_0)$")
            plt.legend()
        elif "keithley" in column:
            assert column in self.data_labels, "column doesn't exist, choose from: {}".format(self.data_labels)
            y,x = self.data[column].index.levels
            I_measure = self.data[column].values
            I_measure = I_measure.reshape(len(I_measure))
            L = len(x) * len(y)
            if len(I_measure) < L:
                I_measure = np.pad(I_measure, (0, L-len(I_measure)))
            try:    
                I_measure = I_measure.reshape((len(y), len(x)))
                y = np.repeat(self.y, len(x)).reshape((len(y), len(x)))
            except ValueError:
                I_measure = I_measure.reshape((len(x), len(y)))
                #y = np.repeat(self.y, len(y)).reshape((len(x), len(y)))
                
            if self.line_resist_sub:
                p = plt.plot(self.V_new[line,:], I_measure[line,:], label = self.y_lab + "= {}".format(round(self.y[line], digs)), **plotkwargs)
            else:
                p = plt.plot(self.x[line,:], I_measure[line,:], label = self.y_lab + "= {}".format(round(self.y[line], digs)), **plotkwargs)
            plt.xlabel(self.x_lab)
            plt.ylabel("I")
            plt.legend()
    
    def plot1D(self, column="lockin 1", **plotkwargs):
        """1d plot of a measured value
        Arguments:
            column: string like 'lockin 1' or 'keithley 2' 
            pcolorkwargs: [OPT] any argument you would put in plt.plot()"""
        
        plt.figure()
        title = self.exp_name + " ({})".format(str(self.id))
        plt.title(title)
        if "lockin" in column:
            assert self.G.any(), "Differential Conductance not computed, run diff_cond"
            #y = np.repeat(self.y, self.G.shape[1]).reshape(self.G.shape)
            if self.line_resist_sub:
                p = plt.plot(self.V_new, self.G, **plotkwargs)
            else:
                p = plt.plot(self.x, self.G, **plotkwargs)
            plt.xlabel(self.x_lab)
            plt.ylabel(r"$dI/dV (G_0)$")
        elif "keithley" in column:
            assert column in self.data_labels, "column doesn't exist, choose from: {}".format(self.data_labels)
            I_measure = self.data[column].values
            index = self.data[column].index
            I_measure = I_measure.reshape(len(I_measure))
                
            if self.line_resist_sub:
                p = plt.plot(self.V_new, I_measure, **plotkwargs)
            else:
                self.x = index.values
                self.x_lab = index.name
                self.y = I_measure
                p = plt.plot(self.x, I_measure, **plotkwargs)
            plt.xlabel(self.x_lab)
            plt.ylabel("I")
        
    def subt_line_resist(self, I_measure_col, lineR=8484, index_bias_2d=1):
        """
        Function for subtracting line resistance from AC and DC voltage
        Arguments:
            I_measure_col: the keithley that measures voltage drop over the device. Note that this is 
            lineR: the line resistance in ohms
            index_bias_2d: which sweep axis is Vbias, either 0 or 1
        creates:
            self.V_new: the new 2D Vbias array. It's 2d because of different resistances per linescan
            self.y, self.y_lab: the y parameter and its label. Often field or gate
        
        """
        assert I_measure_col in self.data_labels, "column doesn't exist, choose from: {}".format(self.data_labels)
        assert not self.line_resist_sub, "Line resist already subtracted"
        I_measure = self.data[I_measure_col].values
        I_measure = I_measure.reshape(len(I_measure))
        if type(self.data[I_measure_col].index) != pd.MultiIndex:
            V_bias = self.data[I_measure_col].index.values
            V_bias_resize = V_bias.reshape(len(I_measure))
        else:
            index_y = 1- index_bias_2d
            V_bias = self.data[I_measure_col].index.levels[index_y].values
            self.y_lab = self.data[I_measure_col].index.levels[index_bias_2d].name
            self.y= self.data[I_measure_col].index.levels[index_bias_2d].values
            L = len(self.y) * len(V_bias)
            if len(I_measure) < L:
                I_measure = np.pad(I_measure, (0, L-len(I_measure)))
            I_measure = I_measure.reshape(len(self.y), (len(V_bias)))
            V_bias_resize = V_bias.reshape(1,len(V_bias))    
        
        R = V_bias_resize/I_measure
        self.V_new = V_bias_resize * (R/(R + lineR))
        for k,v in self.lockins_ac.items():
            _ac = self.lockins_ac[k]
            self.lockins_ac[k] = self.lockins_ac[k] * (R/(R + lineR))
            self.lockins_ac[k][self.lockins_ac[k] ==0] = _ac
        self.x = self.V_new
        self.x_lab = r"$V_{bias}$"
        self.z = I_measure
        self.z_lab = I_measure_col
        self.line_resist_sub = True
    
    def rename_columns(self, dc):
        """rename data columns according to a dictionary
        Arguments:
            dc: dict containing old columns as keys, new columns as values"""
        cols = copy.deepcopy(self.data)
        for k,v in dc.items():
            assert k in self.data_labels, "column doesn't exist, choose from: {}".format(self.data_labels)
            _d = self.data[k] 
            del self.data[k]
            self.data[v] = _d
            
        self.data_labels = self.data.keys()
            
    
    def _rename_instruments_init(self):
        self.lockins_ac= {}
        i = 0
        for instrument in self.instruments:
            a = re.search("lockin", instrument)
            if a:
                i += 1
                V_ac = self.dataset.snapshot['station']["instruments"][instrument]['parameters']['amplitude']['value']
                self.lockins_ac["lockin {}".format(i)] = V_ac
                _d = self.dataset.snapshot['station']["instruments"][instrument] 
                self.dataset.snapshot['station']["instruments"]["lockin {}".format(i)]=_d
                
        self.keithleys= 0
        for instrument in self.instruments:
            a = re.search("keithley", instrument)
            if a:
                self.keithleys += 1
        i=1   
        cols = copy.deepcopy(self.data)
        for col in cols:
            a = re.search("keithley", col)
            if a:
                _d = self.data[col] 
                del self.data[col]
                self.data["keithley {}".format(i)] = _d
                i+=1       
        i=1
        for col in cols:
            a = re.search("K", col)
            if a:
                _d = self.data[col] 
                del self.data[col]
                self.data["keithley {}".format(i)] = _d
                i+=1
        i=1
        for col in cols:
            a = re.search("lockin", col)
            if a:
                _d = self.data[col] 
                del self.data[col]
                self.data["lockin {}".format(i)] = _d
                #self.data["dIdV {}".format(i)] = _d
                i+=1
        for col in cols:
            a = re.search("L", col)
            if a:
                _d = self.data[col] 
                del self.data[col]
                self.data["lockin {}".format(i)] = _d
                #self.data["dIdV {}".format(i)] = _d
                i+=1
        self.data_labels = self.data.keys()
        self.instruments = self.dataset.snapshot['station']["instruments"].keys()

            
class Dataset_qcodes_basic():
    """
    A qcodes dataset class that does the absolute minimum: loads DB, extracts data.
    Intended to spare you from copy pasting cells in Jupyter 
    and should not be broken by however weird measurements.
    Properties:
        dataset (load_by_id)
        dF (pandas.DataFrame)
        run_id; exp_name; sample_name; instruments
        setpoint_labels; data_labels
        params: all qcodes Parameters
        Each Parameter gets its own np.array under the same name
    """
    def __init__(self, data_dir, db_filename, run_id):
        qc.config["core"]["db_location"] = os.path.join(data_dir, db_filename)
        self.dataset = load_by_id(run_id)
        self.dF = self.dataset.get_data_as_pandas_dataframe()
        self.run_id = run_id
        self.exp_name = self.dataset.exp_name
        self.sample_name = self.dataset.sample_name
        self.data_labels = list(self.dF.keys())
        self.setpoint_labels = [p.name for p in self.dataset.get_parameters() \
                                if len(p.depends_on) == 0]
        self.params = self.dataset.get_parameters()
        self.instruments = list(self.dataset.snapshot['station']["instruments"].keys())
        for label in self.setpoint_labels + self.data_labels:
            redund_vals = self.dataset.get_parameter_data(label)[label][label]
            shape = (int(len(redund_vals)/len(self.data_labels)), len(self.data_labels))
            self.__dict__[label] = redund_vals[:shape[0]*shape[1]].reshape(shape)[:,0]
        for label in self.data_labels:
            self.__dict__[label] = np.array(self.dF[label][label].values)
            
    def full_dataset_title(self):
        return f"{self.dataset.path_to_db.split('/')[-1]} ID{self.run_id}: {self.exp_name}"


class Dataset_2d_qcodes(Dataset_qcodes_basic, Dataset_2d):
    """
    A qcodes dataset class that handles 2-dimensional sweeps in a grid:
    loads DB, extracts data, 
    finds out x and y axes, reshapes data into matrices.
    Properties:
        dataset (load_by_id)
        dF (pandas.DataFrame)
        run_id; exp_name; sample_name; instruments
        setpoint_labels; data_labels
        params: all qcodes Parameters
        Each Parameter gets its own np.array under the same name
        x; y; xlabel; ylabel; x_param; y_param
        z_list; zlabel_list; zlabel_dict
        Plotting object dicts: figs; axs; meshes; cbars
    """
    def __init__(self, data_dir, db_filename, run_id):
        Dataset_qcodes_basic.__init__(self, data_dir, db_filename, run_id)
        if len(self.setpoint_labels) != 2:
            print('WARNING: number of independent parameters is not 2!')
        
        def _look_for_param(param_name):
            for p in self.params:
                if p.name == param_name:
                    return p
                
        y_init, x_init = [self.__dict__[label] for label in self.setpoint_labels]
        self.y_param, self.x_param = [_look_for_param(label) \
                                      for label in self.setpoint_labels]
        if y_init[0] != y_init[1]: # makes sure x is the inner/fast axis and y outer/slow
            y_init, x_init = x_init, y_init
            self.y_param, self.x_param = self.x_param, self.y_param
            
        z_list_init = [self.__dict__[label] for label in self.data_labels]
        Dataset_2d.__init__(self, x_init, y_init, z_list_init)
        for kk, label in enumerate(self.data_labels):
            self.__dict__[label] = self.z_list[kk]
        
        self.zlabel_list = []
        self.zlabel_dict = {}
        for label in self.data_labels:
            z_param = _look_for_param(label)
            zlabel = f'{z_param.name} ({z_param.unit})'
            self.zlabel_list.append(zlabel)
            self.zlabel_dict[label] = zlabel
        self.xlabel = f'{self.x_param.name} ({self.x_param.unit})'
        self.ylabel = f'{self.y_param.name} ({self.y_param.unit})'
        
        self.figs, self.axs, self.meshes, self.cbars = {}, {}, {}, {}
        
    def _store_plotting_objects(self, param_name, fig, ax, mesh, cbar):
        self.figs[param_name] = fig
        self.axs[param_name] = ax
        self.meshes[param_name] = mesh
        self.cbars[param_name] = cbar
        
    def imshow_param(self, param_name):
        """
        plt.pcolormesh a measured parameter, without x and y
        Useful when plot_by_id cannot handle non-standard indices
        """
        fig, ax = plt.subplots()
        mesh = ax.imshow(self.__dict__[param_name], cmap='viridis')
        cbar = plt.colorbar(mesh, ax=ax, label=self.zlabel_dict[param_name])
        self._store_plotting_objects(param_name, fig, ax, mesh, cbar)
        
    def plot_param(self, param_name, fig=None, ax=None, cmap='viridis', \
                   vmin=None, vmax=None, gamma=None, sym_cmap=False, \
                   scale_x=1, scale_y=1, scale_z=1):
        """
        Plots the specified parameter name with pcolormesh.
        fig, ax can be either specified or created if not.
        Supports colormap, vmin, vmax, power-law gamma specs.
        If sym_cmap is enabled, the colorbar will be centered at 0 with automatic vmin, vmax.
        Plotting objects are stored in the dicts (figs, axs, meshes, cbars) of the instance
        """
        if fig == None:
            fig, ax = plt.subplots()
        else:
            assert ax != None, 'You supplied a fig to use but no ax'
        
        xlabel = self.xlabel if scale_x == 1 else fr'{scale_x}$\times$'+self.xlabel
        ylabel = self.ylabel if scale_y == 1 else fr'{scale_y}$\times$'+self.ylabel
        try:
            zlabel = self.zlabel_dict[param_name] if scale_z == 1 else f'{scale_z}$\times$'+self.zlabel_dict[param_name]
        except KeyError:
            zlabel = 'No label'
            print('WARNING: You might be plotting a processed value, not measured by qcodes.')
            print('Add a label to it using d.zlabel_dict[var_name] = "var (a.u.)"')
            
        mesh = ax.pcolormesh(self.x*scale_x, self.y*scale_y, \
                             self.__dict__[param_name]*scale_z, cmap=cmap)
        cbar = plt.colorbar(mesh, ax=ax, label=zlabel)
        
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.suptitle(self.full_dataset_title(), fontsize=10)
        
        if sym_cmap:
            vmax = np.abs(self.__dict__[param_name]).max()
            vmin = -vmax
            gamma = 1
        if vmin != None or vmax != None or gamma != None:
            vmin = vmin if vmin != None else self.__dict__[param_name].min()
            vmax = vmax if vmax != None else self.__dict__[param_name].max()
            gamma = gamma if gamma != None else 1
            norm = mpl.colors.PowerNorm(gamma, vmin=vmin, vmax=vmax)
            mesh.set_norm(norm)
            cbar.update_normal(mesh)
        
        self._store_plotting_objects(param_name, fig, ax, mesh, cbar)
        
    def update_minmax(self, param_name, vmin, vmax, gamma=None):
        """
        Looks for the corresponding plotting objects stored in the instance
        and updates the corresponding vmin, vmax, gamma/
        """
        if gamma == None:
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        else:
            # assert gamma > 0, 'Check input: gamma needs to be a positive number'
            norm = mpl.colors.PowerNorm(gamma, vmin=vmin, vmax=vmax)
        self.meshes[param_name].set_norm(norm)
        self.cbars[param_name].update_normal(self.meshes[param_name])


class Dataset_3T_qcodes(Dataset_3T, Dataset_2d_qcodes):
    """
    A 3-terminal data class that interfaces with qcodes DB
    Assumption of measurement: fast axis is V_bias twice, first on the left and then on the right
    Slow axis is some arbitrary parameter.
    Voltage in V, current in A, conductance in G0.
    """
    def __init__(self, data_dir, db_filename, run_id, gain=1e7, fake_sim_meas = True):
        """
        Can only handle two sides having the same gain for now.
        fake_sim_meas option takes the (averaged) center column of the measured half and use it as 
        an estimation of the unmeasured half (when the corresponding lockin excitation is actually off).

        """
        Dataset_2d_qcodes.__init__(self, data_dir, db_filename, run_id)
        self.IL, self.IR = self.__dict__['Current_L']*1e-9, self.__dict__['Current_L']*1e-9
        self.zlabel_dict['IL'] = '$I_L$ (A)'
        self.zlabel_dict['IR'] = '$I_R$ (A)'
        
        self.IL_VL, self.IL_VR = self._split_gXX_in_halves(self.IL)
        self.IR_VL, self.IR_VR = self._split_gXX_in_halves(self.IR)
        self.zlabel_dict['IL_VL'] = '$I_L$ (A)'
        self.zlabel_dict['IL_VR'] = '$I_L$ (A)'
        self.zlabel_dict['IR_VL'] = '$I_R$ (A)'
        self.zlabel_dict['IR_VR'] = '$I_R$ (A)'
        
        for subscript in ['LL', 'LR', 'RL', 'RR']:
            self.__dict__[f'g{subscript}'] = self.__dict__[f'Conductance_{subscript}']
            self.zlabel_dict[f'g{subscript}'] = f'$g_{{{subscript}}}$ ($2e^2/h$)'
            
            self.__dict__[f'g{subscript}_VL'], self.__dict__[f'g{subscript}_VR'] = \
            self._split_gXX_in_halves(self.__dict__[f'g{subscript}'])
            self.zlabel_dict[f'g{subscript}_VL'] = f'$g_{{{subscript}}}$ ($2e^2/h$)'
            self.zlabel_dict[f'g{subscript}_VR'] = f'$g_{{{subscript}}}$ ($2e^2/h$)'
            
        if fake_sim_meas:
            self.gLL_VR = self._fake_useless_half(self.gLL_VL)
            self.gLR_VL = self._fake_useless_half(self.gLR_VR)
            self.gRL_VR = self._fake_useless_half(self.gRL_VL)
            self.gRR_VL = self._fake_useless_half(self.gRR_VR)
        
        Is_VL = [self.IL_VL, self.IR_VL]
        Is_VR = [self.IL_VR, self.IR_VR]
        gs_VL = [self.gLL_VL, self.gLR_VL, self.gRL_VL, self.gRR_VL]
        gs_VR = [self.gLL_VR, self.gLR_VR, self.gRL_VR, self.gRR_VR]
        
        ncols = int(len(self.x)/2)
        if len(self.x) != ncols*2:
            print('WARNING: length of bias vector is an odd number.')
        Dataset_3T.__init__(self, self.x[:ncols], self.x[ncols:], Is_VL, Is_VR, gs_VL, gs_VR, self.y, self.ylabel)
        
    def _split_gXX_in_halves(self, gXX):
        ncols = int(gXX.shape[1] / 2)
        return gXX[:, :ncols], gXX[:, ncols:]
    
    def _fake_useless_half(self, gXX, avg_width = 5):
        ncols = gXX.shape[1]
        gXX_middle_mean = gXX[:, int(ncols*0.5-avg_width*0.5):int(ncols*0.5+avg_width*0.5)].mean(axis=1)
        return gXX_middle_mean[:,np.newaxis].dot(np.ones((1,ncols)))
    
    def update_minmax(self, label, vmin, vmax, gamma=None):
        """
        label can either be in {'LL', 'LR', 'RL', 'RR'} or a param_name
        It calls either the update_minmax in Dataset_3T or Dataset_2d_qcodes accordingly
        """
        if label in {'LL', 'LR', 'RL', 'RR'}:
            Dataset_3T.update_minmax(self, label, vmin, vmax, gamma)
        elif label in self.__dict__:
            Dataset_2d_qcodes.update_minmax(self, label, vmin, vmax, gamma)
        else:
            print('label needs to be either a conductance matrix term such as LL or the name of a property variable')