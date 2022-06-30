

from engens.core.EnGens import EnGen
from engens.core.ClustEn import *
from engens.core.FeatureSelector import *
from engens.core.DimReduction import *
import os

from plotly.offline import iplot
from plotly.graph_objs import graph_objs as go

from ipywidgets import interact, interactive, fixed, interact_manual
from ipywidgets import interactive, fixed, IntSlider, VBox, HBox, Layout
import ipywidgets as widgets
from IPython.display import display,clear_output
from ipywidgets.embed import embed_minimal_html
import matplotlib.colors as mc
import nglview as nv
import time
import pytraj as pt



class PlotUtils(object):

    def __init__(self, engen:EnGen, clust: ClustEn, file_loc: str = "./", stride: int = 1) -> None:
        self.engen = engen
        self.clust = clust
        self.root_loc = file_loc
        self.stride = stride
        if not os.path.exists(self.root_loc):
            os.makedirs(self.root_loc)

        # (x,y) coordinates of the transformed data
        self.x = engen.dimred_data[::self.stride,0]
        self.y =  engen.dimred_data[::self.stride,1]

        # number of cluster representative
        self.rep_fnum = np.array(clust.chosen_frames).flatten()
        
        # (x,y) coordinates of cluster representative
        self.rep_x = engen.dimred_data[::,0][self.rep_fnum]
        self.rep_y = engen.dimred_data[::,1][self.rep_fnum]

    def view_timeline(self):
        if self.engen.crystal_flag == True:
            # plot for crystals
            pfig = go.Figure()
            clust_data = self.clust.labels[self.clust.chosen_index]
            clust_data_chosen = clust_data[::self.stride]

            x_labs = [x for x in self.engen.pdb_files]

            scatter = go.Scatter(x=x_labs, y=clust_data_chosen,              
                    marker=dict(
                        size=5,
                        color=clust_data_chosen,
                        colorscale="Viridis"
                    ),
                mode="markers"
                )
            pfig.add_trace(scatter)
            colors = cm.viridis(np.linspace(0, self.clust.chosen_index+2, self.clust.chosen_index+2)/(self.clust.chosen_index+2))
            for i, rf in enumerate(self.rep_fnum): 
                atxt = "<b>C"+str(self.clust.labels[self.clust.chosen_index][rf])+"; F"+str(rf)+"</b>"
                pfig.add_vline(x=rf, line_width=1, line_dash="dash", line_color="red")
                pfig.add_annotation(text=atxt,
                        x=rf, y=1.06, showarrow=False,
                        xref="x", yref="paper",
                        font=dict(
                                color="rgb"+str(mc.to_rgb(colors[self.clust.labels[self.clust.chosen_index][rf]][:3]))
                        ),
                        textangle=20)
            pfig.update_layout(go.Layout(width=800, height=400),
            xaxis_title="pdb_files",
            yaxis_title="cluster",
            hovermode="x"
                        )
            pfig.update_yaxes(dtick="d")
            iplot(pfig)
        else:
            #plot for trajectory
            pfig = go.Figure()
            clust_data = self.clust.labels[self.clust.chosen_index]
            clust_data_chosen = clust_data[::self.stride]

            x_linspace = np.arange(start = 1, stop = clust_data.shape[0], step = self.stride)
            scatter = go.Scatter(x=x_linspace, y=clust_data_chosen,              
                    marker=dict(
                        size=3,
                        color=clust_data_chosen,
                        colorscale="Viridis"
                    ),
                mode="markers"
                )
            pfig.add_trace(scatter)
            colors = cm.viridis(np.linspace(0, self.clust.chosen_index+2, self.clust.chosen_index+2)/(self.clust.chosen_index+2))
            for i, rf in enumerate(self.rep_fnum): 
                atxt = "<b>C"+str(self.clust.labels[self.clust.chosen_index][rf])+"; F"+str(rf)+"</b>"
                pfig.add_vline(x=rf, line_width=1, line_dash="dash", line_color="red")
                pfig.add_annotation(text=atxt,
                        x=rf, y=1.06, showarrow=False,
                        xref="x", yref="paper",
                        font=dict(
                                color="rgb"+str(mc.to_rgb(colors[self.clust.labels[self.clust.chosen_index][rf]][:3]))
                        ),
                        textangle=20)
            pfig.update_layout(go.Layout(width=800, height=400),
            xaxis_title="trajectory frames",
            yaxis_title="cluster",
            hovermode="x"
                        )
            pfig.update_yaxes(dtick="d")
            iplot(pfig)

    def view_PCs(self):
        marker_labels = []
        clust_data = self.clust.labels[self.clust.chosen_index]
        clust_data_chosen = clust_data if self.stride == 0 else clust_data[::self.stride]
        lablist = list(clust_data_chosen)
        frame_ids = np.arange(0,clust_data.shape[0], step=self.stride)

        for i, elem in enumerate(lablist):
            if self.engen.crystal_flag == True:
                marker_labels.append("cluster="+str(elem)+" \n pdb file ="+str(self.engen.pdb_files[i]))
            else:
                marker_labels.append("cluster="+str(elem)+" \n frame="+str(frame_ids[i]))
        
        pfig = go.Figure()
        scatter = go.Scatter(x=self.x, y=self.y,              
                marker=dict(
                    size=3,
                    color=clust_data_chosen,
                    colorscale="Viridis"
                ),
            mode="markers",
            text = marker_labels,
            hovertemplate = "<b>%{text}</b>"
            )
        pfig.add_trace(scatter)
        scatter = go.Scatter(x=self.rep_x, y=self.rep_y,              
                marker=dict(
                    size=6,
                    color='red',
                    line=dict(width=2,color='DarkSlateGrey')
                ),
            mode="markers",
            text = self.rep_fnum,
            hovertemplate = "<b>Representative (frame - %{text})</b>"
            )
        pfig.add_trace(scatter)
        pfig.update_layout(go.Layout(width=450, height=450,showlegend=False),
        xaxis_title="C1",
        yaxis_title="C2")
        iplot(pfig)

    def view_bar(self):
        bar = go.Bar(x=list(range(self.clust.chosen_index+2)) ,
                    y=self.clust.cluster_weights(self.clust.chosen_index),
                    marker={'color': list(range(self.clust.chosen_index+2)), 'colorscale': "Viridis"})
        pfig = go.Figure()
        pfig.add_trace(bar)
        pfig.update_layout(go.Layout(width=450, height=450),
            xaxis_title="cluster",
            yaxis_title="cluster weight")
        pfig.update_xaxes(dtick="d")
        if not self.clust.thr == None:
            pfig.add_hline(y=self.clust.thr, line_width=2, line_dash="dash", line_color="black", 
                    annotation_text="cluster cutoff")
        iplot(pfig)

    def view_extracted_structures(self, path):
        colors = cm.viridis(np.linspace(0, self.clust.chosen_index+2, self.clust.chosen_index+2)/(self.clust.chosen_index+2))
        
        res_pdbs = []
        for file in os.listdir(path):
            if file[-4:] == ".pdb":
                res_pdbs.append(file)

        if not len(res_pdbs)==len(self.clust.chosen_frames):
            print("There is a mismatch between the number of PDB files in the directory\n you provided and the number of selected conformations.")
            print("Structures will not be displayed.")
            print("Make sure you provide the directory you provide \n contains only EnGens generated conformations.")
            return None

        nw = nv.NGLWidget()
        for i, pdb_file in enumerate(res_pdbs):
            clust_n = self.clust.labels[self.clust.chosen_index][self.clust.chosen_frames[i]] 
            name = "Cluster - "+str(clust_n)
            nw.add_component(os.path.join(path,pdb_file), name=name, default_representation=False, ext="pdb")
            nw[i].add_cartoon(color = mc.rgb2hex(colors[clust_n]))
        nw.center()
        return nw

    def view_movie(self, repr = "cartoon"):
        cluster_colors = cm.viridis(np.linspace(0, self.clust.chosen_index+2, self.clust.chosen_index+2)/(self.clust.chosen_index+2))
        frame_colors = cluster_colors[self.clust.labels[self.clust.chosen_index]]
        def on_frame_change(change):
            frame = change['new']
            c = mc.rgb2hex(frame_colors[frame])
            update_func = getattr(view, "update_"+repr)
            update_func(color = c)
            time.sleep(0.01)
        view = nv.show_pytraj(pt.load(self.engen.full_traj_name, self.engen.full_ref))
        view.clear_representations()
        key = "add_"+repr
        add_func = getattr(view, key)
        add_func()
        update_func = getattr(view, "update_"+repr)
        update_func(color = mc.rgb2hex(frame_colors[0]))

        view.observe(on_frame_change, names=['frame'])
        return view

    def dashboard1(self, path = None):
        widg1 = interactive(self.view_timeline)
        widg2 = interactive(self.view_PCs)
        widg3 = interactive(self.view_bar)
        widg4 = self.view_extracted_structures(path)
        if widg4 == None:
            display(VBox([HBox([widg2, widg3]), widg1]))
        else:
            display(VBox([HBox([widg2, widg3]), widg1, widg4]))
