import contextlib
from engens.core.EnGens import EnGen
import pyemma
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
from hde import HDE
from math import sqrt
import tensorflow as tf


class DimReduction(object):

    def __init__(self) -> None:
        self.transformed_data = None
        self.component_number = None
        self.reducer = None
        self.engen = None
        super().__init__()

    def apply(self):
        #implement feature selection
        self.engen.dimred_data = self.transformed_data
        self.engen.data = self.transformed_data

    def choose_n(self, n:int):
        self.component_number = n
        self.transformed_data = self.reducer.get_output(dimensions=np.arange(n))[0]




class PCAReducer(DimReduction):

    def __init__(self, engen:EnGen) -> None:
        
        super().__init__()
        if engen.data is None:
            raise Exception("No data generated with this EnGen!")
            return
        if engen.chosen_feat_index == -1:
            raise Exception("Features not chosen yet!")
            return
        self.engen = engen
        self.data = self.engen.data[self.engen.chosen_feat_index][1]
        self.reducer = pyemma.coordinates.pca(self.data)
        self.transformed_data = self.reducer.get_output()[0]


    def update_component_number(self, number:int):
        self.component_number = number
        self.transformed_data = self.reducer.get_output(dimensions=np.arange(number))[0]

    
    def plot_2d(self, save_loc:str=None) -> None:
        
        y = self.transformed_data
        pyemma.plots.plot_free_energy(y[:,0], y[:,1], cbar=True)
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        if not save_loc is None: plt.savefig(save_loc)

    def plot_3d(self, save_loc:str=None) -> None:
        # Plot PC1, PC2, PC3
        pca_y = self.transformed_data
        fig = px.scatter_3d(x=pca_y[:,0], y=pca_y[:,1], z=pca_y[:,2])
        fig.update_traces(marker=dict(size=1,
                                    line=dict(width=0,
                                                color='DarkSlateGrey')),
                        selector=dict(mode='markers'))
        fig.show()
        if not save_loc is None: fig.write_image(save_loc)

    def plot_variance(self, var_thr=90, save_loc:str=None)->None:
        pca_eigenvalues = self.reducer.eigenvalues
        variance = np.cumsum(pca_eigenvalues)/np.sum(pca_eigenvalues) * 100
        pca_num = 0
        for i, v in enumerate(variance):
            if v >= var_thr:
                pca_num = i
                break
        plt.scatter(np.arange(len(pca_eigenvalues)), variance)
        plt.axvline(pca_num, color='red')
        plt.xlabel("Principal Components (Eigenvalue Index)")
        plt.ylabel("Variance explained (%) (Eigenvalue)")
        plt.title("PCA explained variance (with thr = {})".format(var_thr))
        if not save_loc is None: plt.savefig(save_loc)

        print("Total of "+str(v)+"% of variance explaned by first "+str(pca_num)+" PCs.")

    def get_variance(self, var_thr=90, save_loc:str=None)->None:
        pca_eigenvalues = self.reducer.eigenvalues
        variance = np.cumsum(pca_eigenvalues)/np.sum(pca_eigenvalues) * 100
        pca_num = 0
        for i, v in enumerate(variance):
            if round(v, 4) >= var_thr:
                pca_num = i
                break
        print("PCA explained variance (with thr = {}) by first {} components".format(var_thr, pca_num))
        return pca_num



class TICAReducer(DimReduction):

    def __init__(self, engen:EnGen, TICA_lagtimes:list=[1,2, 5, 10, 25, 50]) -> None:
        
        super().__init__()
        
        if engen.data is None:
            raise Exception("No data generated with this EnGen!")
            return
        if engen.chosen_feat_index == -1:
            raise Exception("Features not chosen yet!")
            return
        self.engen = engen
        self.data = self.engen.data[self.engen.chosen_feat_index][1]
        self.TICA_lagtimes = TICA_lagtimes
        print("Transforming with TICA - might take some time!")
        self.tica_objs = []
        self.tica_objs_ts = []
        for l in TICA_lagtimes:
            print("lag:",l)
            tica_obj = pyemma.coordinates.tica(self.data, lag=l, var_cutoff=1)
            self.tica_objs.append(tica_obj)
            tica_timescales = tica_obj.timescales
            self.tica_objs_ts.append(tica_timescales)
        self.tica_obj = None
        self.transformed_data = None
        self.component_number = None

    def update_component_number(self, number:int):
        self.component_number = number
        self.transformed_data = self.tica_obj.get_output(dimensions=np.arange(number))[0]

    def plot_lag_analysis(self, timescale_num:int=10, chosen_lag:int=50, save_loc:str=None):
        # Plot tica timescales as a function of lag time
        lags = self.TICA_lagtimes
        nlags = len(lags)
        timescale_num = min(timescale_num, self.tica_objs_ts[0].shape[0])
        ts_list = np.zeros((nlags, timescale_num))
        for i, lag in enumerate(lags):
            timescales = self.tica_objs_ts[i]
            ts_list[i, :] = timescales[:timescale_num]


        plt.semilogy(lags, ts_list)
        plt.ylabel('Timescales (ns)')
        plt.xlabel('Lag time (ns)')
        plt.fill_between(lags, 1, lags, facecolor='Gray')
        plt.axvline(chosen_lag, linewidth=2, color='black')
        if not save_loc is None: plt.savefig(save_loc)

    def choose_lag(self, lag:int, tic_thr:int=None):
        self.tica_obj = pyemma.coordinates.tica(self.data, lag=lag, var_cutoff=1)
        self.reducer = self.tica_obj
        self.transformed_data = self.tica_obj.get_output()[0]

    def resolved_processes(self, timescales, lag):
        tmp_diff = np.sort(np.abs(np.diff(timescales)))[::-1]
        signifficant_thr = np.mean(tmp_diff)+2*np.std(tmp_diff)
        diff = np.abs(np.diff(timescales))
        resolved_p = []
        for i, elem in enumerate(list(diff)):
            if elem < signifficant_thr: break
            if elem < lag: break
            resolved_p.append((i,elem))
        return resolved_p

    def choose_lag_auto(self):
        lag_num = len(self.TICA_lagtimes)
        best_lags = []
        #for all processes 
        for i in range(lag_num):
            process_ts = [ts[i] for ts in self.tica_objs_ts]

            #if process time scale is below lagtime - do not count it
            if process_ts[-1] < self.TICA_lagtimes[-1]:
                break 

            x1, y1 = self.TICA_lagtimes[0], process_ts[0]
            x2, y2 = self.TICA_lagtimes[lag_num-1], process_ts[lag_num-1]
            distances = []
            for j in range(lag_num):
                x0 = self.TICA_lagtimes[j]
                y0 = process_ts[j]
                numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
                denominator = sqrt((y2 - y1)**2 + (x2 - x1)**2)
                distances.append(numerator/denominator)
            optimal_j = distances.index(max(distances))
            optimal_lag = self.TICA_lagtimes[optimal_j]
            best_lags.append(optimal_lag)

        if len(best_lags) == 0:
            print("No lags found appropriate - adjust your list of lags (make upper bound lower).")
            return None

        best_lag = int(min(best_lags))
        best_lag_index = self.TICA_lagtimes.index(best_lag)
        print("Chosen lag time: {}".format(best_lag))
        self.tica_obj = pyemma.coordinates.tica(self.data, lag=best_lag, var_cutoff=1)
        self.reducer = self.tica_obj
        self.transformed_data = self.tica_obj.get_output()[0]
        res_proc_n = i
        print("Number of processes above lag time: {}".format(res_proc_n))
        res_ps = self.resolved_processes(self.tica_objs_ts[best_lag_index], best_lag_index)
        print("Number of clearly resolved processes with TICA: {}".format(len(res_ps)))
        print("Processes (index, timescale): ")
        print(res_ps)
        self.choose_lag(best_lag)
        return best_lag
    
    def plot_2d(self, save_loc:str=None) -> None:
        
        if self.tica_obj is None: raise Exception("Lag not chosen!")
        y = self.transformed_data
        Y_concat = y
        pyemma.plots.plot_free_energy(Y_concat[:,0], Y_concat[:,1], cbar=True)
        plt.xlabel("TIC1")
        plt.ylabel("TIC2")
        if not save_loc is None: plt.savefig(save_loc)

    def plot_3d(self, save_loc:str=None) -> None:
        if self.tica_obj is None: raise Exception("Lag not chosen!")
        # Plot PC1, PC2, PC3
        y = self.transformed_data
        fig = px.scatter_3d(x=y[:,0], y=y[:,1], z=y[:,2], color=list(range(y.shape[0])))
        fig.update_traces(marker=dict(size=1,
                                    line=dict(width=0,
                                                color='DarkSlateGrey')),
                        selector=dict(mode='markers'))
        fig.layout.coloraxis.colorbar.title = 'frame number'
        fig.show()
        if not save_loc is None: fig.write_image(save_loc)

    def plot_variance(self, var_thr=90, save_loc:str=None)->None:
        
        if self.tica_obj is None: raise Exception("Lag not chosen!")
        tica_eigenvalues = self.tica_obj.eigenvalues
        variance = np.cumsum(tica_eigenvalues**2)/ np.sum(tica_eigenvalues**2)*100
        tica_num = 0
        for i, v in enumerate(variance):
            if v >= var_thr:
                tica_num = i
                break
                
        plt.scatter(np.arange(len(tica_eigenvalues)), variance)
        plt.xlabel("TICA Components (Eigenvalue Index)")
        plt.ylabel("Kinetic variance explained (%) (Eigenvalue^2)")
        plt.axvline(tica_num, color='red')
        print("Total of "+str(v)+"% of variance explaned by first "+str(tica_num)+" ICs.")
        if not save_loc is None: plt.savefig(save_loc)
        
    def get_variance(self, var_thr=90, save_loc:str=None)->None:
        
        if self.tica_obj is None: raise Exception("Lag not chosen!")
        tica_eigenvalues = self.tica_obj.eigenvalues
        variance = np.cumsum(tica_eigenvalues**2)/ np.sum(tica_eigenvalues**2)*100
        tica_num = 0
        for i, v in enumerate(variance):
            if round(v, 4) >= var_thr:
                tica_num = i
                break
        print("TICA explained variance (with thr = {}) by first {} components".format(var_thr, tica_num))
        return tica_num
                

class HDEReducer(DimReduction):

    def __init__(self, engen:EnGen, HDE_lagtimes:list=[1,2, 5, 10, 25, 50], n_comp=20) -> None:
        
        super().__init__()
        
        tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

        if engen.data is None:
            raise Exception("No data generated with this EnGen!")

        if engen.chosen_feat_index == -1:
            raise Exception("Features not chosen yet!")

        self.engen = engen
        self.data = self.engen.data[self.engen.chosen_feat_index][1]
        self.HDE_lagtimes = HDE_lagtimes
        self.n_components = n_comp
        print("Transforming with HDE - might take some time!")
        self.hde_objs = []
        self.hde_objs_ts = []
        for l in HDE_lagtimes:
            print("lag:",l)
            print("number of components:",n_comp)
            print("feat shape:", self.data.shape[1])
            print("data shape:", self.data.shape[0])
            model = HDE(
                self.data.shape[1], 
                n_components=n_comp, 
                n_epochs=20, 
                lag_time=l,
                batch_size=self.data.shape[0],
                batch_normalization=False
            )
                
            with open('HDE_log.txt','a') as f:
                with contextlib.redirect_stdout(f):
                    model.fit_transform(self.data)
            timescales = model.timescales_
            self.hde_objs.append(model)
            hde_timescales = timescales
            self.hde_objs_ts.append(hde_timescales)
        self.hde_obj = None
        self.transformed_data = None

    def plot_lag_analysis(self, chosen_lag:int=50, save_loc:str=None):
        # Plot tica timescales as a function of lag time
        lags = self.HDE_lagtimes
        nlags = len(lags)
        ts_list = np.zeros((nlags, self.n_components))
        for i, lag in enumerate(lags):
            timescales = self.hde_objs_ts[i]
            ts_list[i, :] = timescales[:self.n_components]


        plt.semilogy(lags, ts_list)
        plt.ylabel('Timescales (ns)')
        plt.xlabel('Lag time (ns)')
        plt.fill_between(lags, 1, lags, facecolor='Gray')
        plt.axvline(chosen_lag, linewidth=2, color='black')
        if not save_loc is None: plt.savefig(save_loc)

    def choose_lag(self, lag:int, n_comp:int=20):
        self.hde_obj = HDE(
                self.data.shape[1], 
                n_components=n_comp, 
                n_epochs=20, 
                lag_time=lag,
                batch_size=self.data.shape[0],
                batch_normalization=False
            )
        self.reducer = self.hde_obj
        with open('HDE_log.txt','a') as f:
            with contextlib.redirect_stdout(f):
                self.transformed_data = self.hde_obj.fit_transform(self.data)
    
    
    def resolved_processes(self, timescales, lag):
        tmp_diff = np.sort(np.abs(np.diff(timescales)))[::-1]
        signifficant_thr = np.mean(tmp_diff)+2*np.std(tmp_diff)
        diff = np.abs(np.diff(timescales))
        resolved_p = []
        for i, elem in enumerate(list(diff)):
            if elem < signifficant_thr: break
            if elem < lag: break
            resolved_p.append((i,elem))
        return resolved_p

    def choose_lag_auto(self):
        best_lags = []
        for i in range(self.n_components):
            process_ts = np.array(self.hde_objs_ts)[:,i]
            lag_num = len(process_ts[~np.isnan(process_ts)])
            # if the process does not have time scales throughout 
            # the whole lag space it was not resolved and should not be used
            if lag_num < len(self.HDE_lagtimes):
                break
            
            x1, y1 = self.HDE_lagtimes[0], process_ts[0]
            x2, y2 = self.HDE_lagtimes[lag_num-1], process_ts[lag_num-1]
            distances = []
            for j in range(len(process_ts)):
                x0 = self.HDE_lagtimes[j]
                y0 = process_ts[j]
                numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
                denominator = sqrt((y2 - y1)**2 + (x2 - x1)**2)
                distances.append(numerator/denominator)
            optimal_j = distances.index(max(distances))
            optimal_lag = self.HDE_lagtimes[optimal_j]
            best_lags.append(optimal_lag)

        if len(best_lags)==0: 
            print("No processes resolved, no best lag, consider changing HDE output dimensions.")
            return None
        res_proc_num = i
        best_lag = int(min(best_lags))
        best_lag_index = self.HDE_lagtimes.index(best_lag)
        print("Chosen lag time: {}".format(best_lag))
        res_ps = self.resolved_processes(self.hde_objs_ts[best_lag_index], best_lag_index)
    
        print("Number of clearly resolved processes with HDE: {}".format(len(res_ps)))
        print("Processes (index, timescale): ")
        print(res_ps)

        self.hde_obj = HDE(
                self.data.shape[1], 
                n_components= self.n_components, 
                n_epochs=20, 
                lag_time=best_lag,
                batch_size=self.data.shape[0],
                batch_normalization=False
            )
        self.reducer = self.hde_obj
        with open('HDE_log.txt','a') as f:
            with contextlib.redirect_stdout(f):
                self.transformed_data = self.hde_obj.fit_transform(self.data)
        return best_lag
    
    def plot_2d(self, save_loc:str=None) -> None:
        
        if self.hde_obj is None: raise Exception("Lag not chosen!")
        y = self.transformed_data
        pyemma.plots.plot_free_energy(y[:,0], y[:,1], cbar=True)
        plt.xlabel("HDE-C1")
        plt.ylabel("HDE-C2")
        if not save_loc is None: plt.savefig(save_loc)

    def plot_3d(self, save_loc:str=None) -> None:
        if self.hde_obj is None: raise Exception("Lag not chosen!")
        # Plot PC1, PC2, PC3
        y = self.transformed_data
        fig = px.scatter_3d(x=y[:,0], y=y[:,1], z=y[:,2], color=list(range(y.shape[0])))
        fig.update_traces(marker=dict(size=1,
                                    line=dict(width=0,
                                                color='DarkSlateGrey')),
                        selector=dict(mode='markers'))
        fig.layout.coloraxis.colorbar.title = 'frame number'
        fig.show()
        if not save_loc is None: fig.write_image(save_loc)


dimreds = {
    "PCA": PCAReducer,  
    "TICA": TICAReducer,
    "HDE": HDEReducer
}


