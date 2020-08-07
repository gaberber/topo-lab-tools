import os, re
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from qcodes import load_by_run_spec
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


def batch_guess_2D_dims(inner_axis, outer_axis, data_2d_list, rescale_xy = True):
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
  
class Dataset_qcodes(Dataset_2d):
    def __init__(self, run_id):
        """
        Reads a qcodes database.
        """
        self.dataset = load_by_run_spec(captured_run_id=run_id)
        self.id = run_id
        self.exp_name = self.dataset.exp_name
        self.sample_name = self.dataset.sample_name
        self.data = self.dataset.get_data_as_pandas_dataframe()
        self.data_labels = self.data.keys()
        self.instruments = self.dataset.snapshot['station']["instruments"].keys()
        self.rename_instruments()
        self.line_restist_sub = False
        
    def diff_cond(self, lockin=1, ac_scaling = 1):
        #G = self.data["dIdV {}".format(i)]
        if not self.line_restist_sub:
            print("Beware, line resistance is not subtracted yet")
        dI = self.data["lockin {}".format(lockin)]
        dV = ac_scaling * self.lockins_ac["lockin {}".format(lockin)]
        self.G = 12906*dI.values.reshape(len(dI))/dV.reshape(len(dI))
        
                
    def subt_line_resist(self, I_measure_col, lineR=8484):
        assert I_measure_col in self.data_labels, "column doesn't exist, choose from: {}".format(self.data_labels)
        assert not self.line_restist_sub, "Line resist already subtracted"
        I_measure = self.data[I_measure_col].values
        I_measure = I_measure.reshape(len(I_measure))
        V_bias = self.data[I_measure_col].index.values
        V_bias = V_bias.reshape(len(I_measure))
        R = V_bias/I_measure
        self.V_new = V_bias * (R/(R + lineR))
        for k,v in self.lockins_ac.items():
            self.lockins_ac[k] = self.lockins_ac[k] * (R/(R + lineR))
        
        self.line_restist_sub = True
            
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

            
        