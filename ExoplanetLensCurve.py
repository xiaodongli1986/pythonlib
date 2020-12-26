
import numpy as np
import pandas as pd
make_plot = True
if make_plot:
   import matplotlib.pyplot as plt

class boma3:
    def __init__(self, filepath = '/home/xiaodongli/data/exoplanets/'):
        self._binfile = filepath+'boma3.bin'
        self._paramsfile = filepath+'boma3_params.bin'
        self._param_names = ['planet_star_distance', 'planet_star_massrat', 'angle']
        self._DESCRIB = '''
        
        Usage: 
             boma = boma3() ## initialize the class
             #X, y = boma.data  ## X, y 
             trival_curve = boma.trival_curve ## "trival" curve with no secondary peak; 
                 # mean of the 1476 curves classified to "trival" classes via a MiniBatch clustering
             distance_to_trival = boma3.distance_to_trival  # distance between the curve and the trival curve
        
        There are 3 parameters varying
            planet_star_distance: 
                angular distance between planet and star (normalized by Einstein ring)
            planet star mass-rat
                mass rat between planet and star
            angle 
                angle between planet-star connection and motion between the lens and the background star
        '''
        
        self._trival_curve_mean_file = filepath+'trival_curve_mean.txt'
        self._trival_curve_ids_file = filepath+'trival_curve_ids.txt'
        self._distance_to_trival_center_file = filepath+'distance_to_trival_center.txt'
        self.data = [np.fromfile(self._binfile, dtype = float).reshape(10000, 700),
                     np.fromfile(self._paramsfile, dtype = float).reshape(10000, 7)[:,3:6]]
        
        self.trival_curve, self._trival_curve_ids, self.distance_to_trival = \
            np.loadtxt(self._trival_curve_mean_file), np.loadtxt(self._trival_curve_ids_file),\
            np.loadtxt(self._distance_to_trival_center_file)
        #print(self.distance_to_trival.shape)
        self.df = pd.DataFrame(np.column_stack([
            self.distance_to_trival, np.fromfile(self._paramsfile, dtype = float).reshape(10000, 7)[:,3:6],]),
            columns = ['id', 'trival_rank', 'trival_distance',  ] + self._param_names )
        self.df['id'] = self.df['id'].map(int); 
        self.df['trival_rank'] = self.df['trival_rank'].map(int)
        self.df = self.df.sort_values('trival_rank')
    def rank_select(self, id1, id2, make_plot=False, plot_i1=None, plot_i2=None):
        '''selecting curves whose distance-to-trival distributes between id1-id2
	Usage: X, y = boma.rank_select(7000, 7030, make_plot=True) ## select curve 7000-7030; make a plot'''
        ids = np.array(self.df['id'][id1:id2])
        X, y = self.data[0][ids], self.data[1][ids]
        if make_plot:
            plot_x = np.array(list(range(len(X[0]))))
            if plot_i1 == None or plot_i2 == None:
               plot_i1, plot_i2 = 0, len(plot_x)
            fig, ax = plt.subplots(figsize=(14,7))
            ax.plot(plot_x[plot_i1:plot_i2], self.trival_curve[plot_i1:plot_i2], lw=3, ls='--', c='k')
            for xx in X:
                ax.plot(plot_x[plot_i1:plot_i2], xx[plot_i1:plot_i2], lw=0.5)
            ax.set_title('curves '+str(id1)+' to '+str(id2))
        return X, y
    
