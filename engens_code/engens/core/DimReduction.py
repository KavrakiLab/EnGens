from engens.core.EnGens import EnGen
import pyemma
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
from hde import HDE


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

    def choose_n(self, n:int):
        self.component_number = n
        self.transformed_data = self.reducer.get_output(dimensions=np.arange(n))[0]




class PCAReducer(DimReduction):

    def __init__(self, engen:EnGen) -> None:
        
        super().__init__()
        if engen.data == None:
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
        if not save_loc == None: plt.savefig(save_loc)

    def plot_3d(self, save_loc:str=None) -> None:
        # Plot PC1, PC2, PC3
        pca_y = self.transformed_data
        fig = px.scatter_3d(x=pca_y[:,0], y=pca_y[:,1], z=pca_y[:,2])
        fig.update_traces(marker=dict(size=1,
                                    line=dict(width=0,
                                                color='DarkSlateGrey')),
                        selector=dict(mode='markers'))
        fig.show()
        if not save_loc == None: fig.write_image(save_loc)

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
        if not save_loc == None: plt.savefig(save_loc)

        print("Total of "+str(v)+"% of variance explaned by first "+str(pca_num)+" PCs.")



class TICAReducer(DimReduction):

    def __init__(self, engen:EnGen, TICA_lagtimes:list=[1,2, 5, 10, 25, 50]) -> None:
        
        super().__init__()
        
        if engen.data == None:
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
            tica_obj = pyemma.coordinates.tica(self.data, lag=l)
            self.tica_objs.append(tica_obj)
            tica_timescales = tica_obj.timescales
            self.tica_objs_ts.append(tica_timescales)
        self.tice_obj = None
        self.transformed_data = None
        self.component_number = None

    def update_component_number(self, number:int):
        self.component_number = number
        self.transformed_data = self.tica_obj.get_output(dimensions=np.arange(number))[0]

    def plot_lag_analysis(self, timescale_num:int=10, chosen_lag:int=50, save_loc:str=None):
        # Plot tica timescales as a function of lag time
        lags = self.TICA_lagtimes
        nlags = len(lags)
        ts_list = np.zeros((nlags, timescale_num))
        for i, lag in enumerate(lags):
            timescales = self.tica_objs_ts[i]
            ts_list[i, :] = timescales[:timescale_num]


        plt.semilogy(lags, ts_list)
        plt.ylabel('Timescales (ns)')
        plt.xlabel('Lag time (ns)')
        plt.fill_between(lags, 1, lags, facecolor='Gray')
        plt.axvline(chosen_lag, linewidth=2, color='black')
        if not save_loc == None: plt.savefig(save_loc)

    def choose_lag(self, lag:int, tic_thr:int=None):
        self.tica_obj = pyemma.coordinates.tica(self.data, lag=lag)
        self.reducer = self.tica_obj
        self.transformed_data = self.tica_obj.get_output()[0]

    
    def plot_2d(self, save_loc:str=None) -> None:
        
        if self.tica_obj == None: raise Exception("Lag not chosen!")
        y = self.transformed_data
        Y_concat = y
        pyemma.plots.plot_free_energy(Y_concat[:,0], Y_concat[:,1], cbar=True)
        plt.xlabel("TIC1")
        plt.ylabel("TIC2")
        if not save_loc == None: plt.savefig(save_loc)

    def plot_3d(self, save_loc:str=None) -> None:
        if self.tica_obj == None: raise Exception("Lag not chosen!")
        # Plot PC1, PC2, PC3
        y = self.transformed_data
        fig = px.scatter_3d(x=y[:,0], y=y[:,1], z=y[:,2], color=list(range(y.shape[0])))
        fig.update_traces(marker=dict(size=1,
                                    line=dict(width=0,
                                                color='DarkSlateGrey')),
                        selector=dict(mode='markers'))
        fig.layout.coloraxis.colorbar.title = 'frame number'
        fig.show()
        if not save_loc == None: fig.write_image(save_loc)

    def plot_variance(self, var_thr=90, save_loc:str=None)->None:
        
        if self.tica_obj == None: raise Exception("Lag not chosen!")
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
        if not save_loc == None: plt.savefig(save_loc)
                

class HDEReducer(DimReduction):

    def __init__(self, engen:EnGen, HDE_lagtimes:list=[1,2, 5, 10, 25, 50], n_comp=20) -> None:
        
        super().__init__()
        
        if engen.data == None:
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

            model = HDE(
                self.data.shape[1], 
                n_components=n_comp, 
                n_epochs=20, 
                lag_time=l,
                batch_normalization=True
            )

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
        if not save_loc == None: plt.savefig(save_loc)

    def choose_lag(self, lag:int, n_comp:int=20):
        self.hde_obj = HDE(
                self.data.shape[1], 
                n_components=n_comp, 
                n_epochs=20, 
                lag_time=lag,
                batch_normalization=True
            )
        self.reducer = self.hde_obj
        self.transformed_data = self.hde_obj.fit_transform(self.data)

    
    def plot_2d(self, save_loc:str=None) -> None:
        
        if self.hde_obj == None: raise Exception("Lag not chosen!")
        y = self.transformed_data
        pyemma.plots.plot_free_energy(y[:,0], y[:,1], cbar=True)
        plt.xlabel("HDE-C1")
        plt.ylabel("HDE-C2")
        if not save_loc == None: plt.savefig(save_loc)

    def plot_3d(self, save_loc:str=None) -> None:
        if self.hde_obj == None: raise Exception("Lag not chosen!")
        # Plot PC1, PC2, PC3
        y = self.transformed_data
        fig = px.scatter_3d(x=y[:,0], y=y[:,1], z=y[:,2], color=list(range(y.shape[0])))
        fig.update_traces(marker=dict(size=1,
                                    line=dict(width=0,
                                                color='DarkSlateGrey')),
                        selector=dict(mode='markers'))
        fig.layout.coloraxis.colorbar.title = 'frame number'
        fig.show()
        if not save_loc == None: fig.write_image(save_loc)


dimreds = {
    "PCA": PCAReducer,  
    "TICA": TICAReducer,
    "HDE": HDEReducer
}


