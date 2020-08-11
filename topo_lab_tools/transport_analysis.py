import os, re
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from qcodes import load_by_run_spec, initialise_or_create_database_at
import qcodes as qc
import pandas as pd
import matplotlib
matplotlib.rcParams["text.usetex"] = False
#matplotlib.rcParams["font.family"] = "Helvetica"
matplotlib.rcParams["text.latex.preamble"] = "\\usepackage{amsmath}"
matplotlib.rcParams["axes.labelsize"] = 16
matplotlib.rcParams["xtick.labelsize"] = 16
matplotlib.rcParams["ytick.labelsize"] = 16
matplotlib.rcParams["font.size"] = 15
matplotlib.rcParams["figure.figsize"] = [12,5]
matplotlib.rcParams["ytick.major.size"] = 8
matplotlib.rcParams["ytick.major.width"] = 1.5
matplotlib.rcParams["ytick.minor.size"] = 4
matplotlib.rcParams["ytick.minor.width"] = 1
matplotlib.rcParams["ytick.direction"] = "in"
matplotlib.rcParams["xtick.major.size"] = 8
matplotlib.rcParams["xtick.major.width"] = 1.5
matplotlib.rcParams["xtick.direction"] = "in"
matplotlib.rcParams["xtick.minor.size"] = 4
matplotlib.rcParams["xtick.minor.width"] = 1
matplotlib.rcParams["savefig.bbox"] = "tight"
matplotlib.rcParams["savefig.pad_inches"] = 0.05
matplotlib.rcParams["axes.linewidth"] = 1.5
matplotlib.rcParams["legend.fontsize"] = 14
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
                print(a[0])
    assert i>0, "'{0}' not found in '{1}' \n These exist: {2}".format(db_name, data_dir, ls)
    initialise_or_create_database_at("{0}/{1}".format(data_dir, db_name) + ".db" )

class Dataset_qcodes(Dataset_2d):
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
            linecut
        """
        
        self.dataset = load_by_run_spec(captured_run_id=run_id)
        self.id = run_id
        self.exp_name = self.dataset.exp_name
        self.sample_name = self.dataset.sample_name
        self.data = self.dataset.get_data_as_pandas_dataframe()
        self.data_labels = self.data.keys()
        self.instruments = self.dataset.snapshot['station']["instruments"].keys()
        self.rename_instruments()
        self.line_resist_sub = False
        
    def diff_cond(self, column, ac_scaling = 1):
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
        fig, ax = plt.subplots(figsize=(15,6))
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
            cbar = plt.colorbar(p)
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
            cbar = plt.colorbar(p)
            cbar.set_label(column)
            
            
    def linecut(self, column="lockin 1", line=0, digs = 2, **plotkwargs):
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
            I_measure = I_measure.reshape(len(I_measure))
                
            if self.line_resist_sub:
                p = plt.plot(self.V_new, I_measure, **plotkwargs)
            else:
                p = plt.plot(self.x, I_measure, **plotkwargs)
            plt.xlabel(self.x_lab)
            plt.ylabel("I")
        
    def subt_line_resist(self, I_measure_col, lineR=8484, index_bias_2d=1):
        """
        Function for subtracting line resistance from AC and DC voltage
        Arguments:
            I_measure_col: the keithley that measures voltage drop over the device. Note that this is 
        
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
            
    def rename_instruments(self):
        self.lockins_ac= {}
        i = 0
        for instrument in self.instruments:
            a = re.search("lockin", instrument)
            if a:
                i += 1
                V_ac = self.dataset.snapshot['station']["instruments"][instrument]['parameters']['amplitude']['value']
                self.lockins_ac["lockin {}".format(i)] = V_ac
                _d = self.dataset.snapshot['station']["instruments"][instrument] 
                #del self.dataset.snapshot['station']["instruments"][instrument]
                #print(self.dataset.snapshot['station']["instruments"].keys())
                self.dataset.snapshot['station']["instruments"]["lockin {}".format(i)]=_d
                
        self.keithleys= 0
        for instrument in self.instruments:
            a = re.search("keithley", instrument)
            if a:
                self.keithleys += 1
        i=1   
        cols = self.data.keys()
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

            
        