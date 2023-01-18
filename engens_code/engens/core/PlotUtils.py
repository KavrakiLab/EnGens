

from engens.core.EnGens import EnGen
from engens.core.ClustEn import *
from engens.core.FeatureSelector import *
from engens.core.DimReduction import *
import os

from plotly.offline import iplot
from plotly.graph_objs import graph_objs as go

from ipywidgets import VBox, HBox
import ipywidgets as widgets
from IPython.display import display
import matplotlib.colors as mc
import nglview as nv
import time
import pytraj as pt
import mdtraj
from plotly.express.colors import sample_colorscale


class PlotUtils(object):

    def __init__(self, engen:EnGen, 
                        clust: ClustEn,
                        file_loc: str = "./", 
                        stride: int = 1,
                        colorscale: str = "Viridis") -> None:
        self.engen = engen
        self.clust = clust
        self.root_loc = file_loc
        self.stride = stride
        self.colorscale = colorscale
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

    def view_timeline(self, fileloc=None, width=None) -> go.Figure: 
        if self.engen.crystal_flag == True:
            # plot for crystals
            pfig = go.Figure()

            clust_data = self.clust.labels[self.clust.chosen_index]
            clust_data_chosen = clust_data[::self.stride]

            x_labs = [x for x in self.engen.structure_names]

            scatter = go.Scatter(x=x_labs, y=clust_data_chosen,              
                    marker=dict(
                        size=5,
                        color=clust_data_chosen,
                        colorscale=self.colorscale,
                    line=dict(width=0.5,color='DarkSlateGrey')
                    ),
                mode="markers"
                )
            pfig.add_trace(scatter)
            colors = sample_colorscale(self.colorscale, np.linspace(0, 1, self.clust.chosen_index+2))
            
            for i, rf in enumerate(self.rep_fnum): 
                atxt = "<b>C"+str(self.clust.labels[self.clust.chosen_index][rf])+"; "+str(x_labs[rf])+"</b>"
                pfig.add_vline(x=rf, line_width=1, line_dash="dash", line_color="red")
                pfig.add_annotation(text=atxt,
                        x=rf, y=1.06, showarrow=False,
                        xref="x", yref="paper",
                        font=dict(
                                size=8,
                                color="black"
                        ),
                        textangle=20)

            if width is None: width=800
            pfig.update_layout(go.Layout(width=width, height=400),
            xaxis_title="PDB codes",
            yaxis_title="cluster",
            hovermode="x"
                        )
            pfig.update_yaxes(dtick="d")
            #iplot(pfig)
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
                        colorscale=self.colorscale,
                    line=dict(width=0.5,color='DarkSlateGrey')
                    ),
                mode="markers"
                )
            pfig.add_trace(scatter)
            
            colors = sample_colorscale(self.colorscale, np.linspace(0, 1, self.clust.chosen_index+2))
            for i, rf in enumerate(self.rep_fnum): 
                atxt = "<b>C"+str(self.clust.labels[self.clust.chosen_index][rf])+"; F"+str(rf)+"</b>"
                pfig.add_vline(x=rf, line_width=1, line_dash="dash", line_color="red")
                pfig.add_annotation(text=atxt,
                        x=rf, y=1.06, showarrow=False,
                        xref="x", yref="paper",
                        font=dict(
                                size=8,
                                color="black"
                        ),
                        textangle=20)
            
            if width is None: width=800
            pfig.update_layout(go.Layout(width=width, height=400),
            xaxis_title="trajectory frames",
            yaxis_title="cluster",
            hovermode="x"
                        )
            pfig.update_yaxes(dtick="d")
            #iplot(pfig)
        if not fileloc is None: 
            pfig.write_html(fileloc)
        return pfig

    def view_PCs(self, fileloc=None) -> go.Figure:
        marker_labels = []
        clust_data = self.clust.labels[self.clust.chosen_index]
        clust_data_chosen = clust_data if self.stride == 0 else clust_data[::self.stride]
        lablist = list(clust_data_chosen)
        frame_ids = np.arange(0,clust_data.shape[0], step=self.stride)

        for i, elem in enumerate(lablist):
            if self.engen.crystal_flag == True:
                marker_labels.append("cluster="+str(elem)+" \n pdb file ="+str(self.engen.structure_names[i]))
            else:
                marker_labels.append("cluster="+str(elem)+" \n frame="+str(frame_ids[i]))
        
        rep_labels = []
        rep_colors = []
        colors = sample_colorscale(self.colorscale, np.linspace(0, 1, self.clust.chosen_index+2), colortype="rgb")
        for i, elem in enumerate(self.rep_fnum):
            rep_colors.append(colors[i])
            if self.engen.crystal_flag == True:
                rep_labels.append("cluster="+str(i)+" \n pdb file ="+str(self.engen.structure_names[i]))
            else:
                rep_labels.append("cluster="+str(i)+" \n frame="+str(elem))
        
        pfig = go.Figure()


        scatter = go.Scatter(x=self.x, y=self.y,              
                marker=dict(
                    size=5,
                    color=clust_data_chosen,
                    colorscale=self.colorscale,
                    line=dict(width=0.5,color='DarkSlateGrey')
                ),
            mode="markers",
            text = marker_labels,
            hovertemplate = "<b>%{text}</b>"
            )
        pfig.add_trace(scatter)



        scatter = go.Scatter(x=self.rep_x, y=self.rep_y,              
                marker=dict(
                    size=6,
                    color=rep_colors,
                    line=dict(width=2,color='red')
                ),
            mode="markers",
            text = rep_labels,
            hovertemplate = "<b>Representative (%{text})</b>"
            )
        pfig.add_trace(scatter)
        pfig.update_layout(go.Layout(width=450, height=450,showlegend=False),
        xaxis_title="C1",
        yaxis_title="C2")
        #iplot(pfig)
        if not fileloc is None: pfig.write_html(fileloc)
        return pfig

    def view_bar(self, fileloc=None):
        bar = go.Bar(x=list(range(self.clust.chosen_index+2)) ,
                    y=self.clust.cluster_weights(self.clust.chosen_index),
                    marker={'color': list(range(self.clust.chosen_index+2)), 'colorscale': self.colorscale})
        pfig = go.Figure()
        pfig.add_trace(bar)
        pfig.update_layout(go.Layout(width=450, height=450),
            xaxis_title="cluster",
            yaxis_title="cluster weight")
        pfig.update_xaxes(dtick="d")
        if not self.clust.thr == None:
            pfig.add_hline(y=self.clust.thr, line_width=2, line_dash="dash", line_color="black", 
                    annotation_text="cluster cutoff")
        #iplot(pfig)
        if not fileloc is None: pfig.write_html(fileloc)
        return pfig

    def view_extracted_structures(self, path):
        
        colors = sample_colorscale(self.colorscale, 
        np.linspace(0, 1, self.clust.chosen_index+2), colortype="rgb")
        res_pdbs = []
        for file in os.listdir(path):
            if file[-4:] == ".pdb":
                res_pdbs.append(file)
        res_pdbs = sorted(res_pdbs)
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
            nw[i].add_cartoon(color = colors[clust_n])
        nw.center()
        return nw

    def view_movie(self, repr = "cartoon"):
        cluster_colors = sample_colorscale(self.colorscale, 
        np.linspace(0, 1, self.clust.chosen_index+2), colortype="rgb")
        frame_colors = [cluster_colors[i] for i in self.clust.labels[self.clust.chosen_index]]
        def on_frame_change(change):
            frame = change['new']
            c = frame_colors[frame]
            update_func = getattr(view, "update_"+repr)
            update_func(color = c)
            time.sleep(0.01)
        view = nv.show_pytraj(pt.load(self.engen.full_traj_name, self.engen.full_ref))
        view.clear_representations()
        key = "add_"+repr
        add_func = getattr(view, key)
        add_func()
        update_func = getattr(view, "update_"+repr)
        update_func(color = frame_colors[0])

        view.observe(on_frame_change, names=['frame'])
        return view

    def dashboard1(self, path = None, fileloc=None):
        widg1 = go.FigureWidget(self.view_timeline(fileloc))
        widg2 = go.FigureWidget(self.view_PCs(fileloc))
        widg3 = go.FigureWidget(self.view_bar(fileloc))
        widg4 = self.view_extracted_structures(path)
        if not fileloc == None:
            if widg4 == None:
                dashboard = VBox([HBox([widg2, widg3]), widg1])
            else:
                dashboard = VBox([HBox([widg2, widg3]), widg1, widg4])
            display(dashboard)
        else:
            if widg4 == None:
                display(VBox([HBox([widg2, widg3]), widg1]))
            else:
                display(VBox([HBox([widg2, widg3]), widg1, widg4]))

    def plot_custom_feature_scatter(self, feature, y_title="feature"):
        
        n_clusters = self.clust.chosen_index + 2
        
        colors = sample_colorscale(self.colorscale, np.linspace(0, 1, self.clust.chosen_index+2))

        y_list = []
        x_list = []
        for i in range(n_clusters):
            y = feature[np.argwhere(self.clust.labels[self.clust.chosen_index]==i).flatten()]
            y_list.append(y)
            x_list.append(np.argwhere(self.clust.labels[self.clust.chosen_index]==i).flatten())
        
        print(x_list)
        print(y_list)
        fig = go.Figure()
        for i, y in enumerate(y_list):
            fig.add_trace(go.Scatter(x = x_list[i], y=y, 
            mode='markers', marker={"color":colors[i]}, 
                                                                name="Cluster #{}".format(i)))
        
        clust_data = self.clust.labels[self.clust.chosen_index]
        #frame_ids = np.arange(0,clust_data.shape[0], step=self.stride)
        rep_labels = []
        rep_colors = []
        rep_x = []
        rep_y = []
        colors = sample_colorscale(self.colorscale, np.linspace(0, 1, self.clust.chosen_index+2), colortype="rgb")
        for i, elem in enumerate(self.rep_fnum):
            rep_x.append(elem)
            rep_y.append(feature[elem])
            rep_colors.append(colors[i])
            if self.engen.crystal_flag == True:
                rep_labels.append("cluster="+str(i)+" \n pdb file ="+str(self.engen.structure_names[i]))
            else:
                rep_labels.append("cluster="+str(i)+" \n frame="+str(elem))
        
        scatter = go.Scatter(x=rep_x, y=rep_y,              
                marker=dict(
                    size=6,
                    color=rep_colors,
                    line=dict(width=2,color='red')
                ),
            mode="markers",
            text = rep_labels,
            name = "Representatives",
            hovertemplate = "<b>Representative (%{text})</b>"
            )
        fig.add_trace(scatter)

        fig.update_layout(yaxis_title=y_title)
        fig.show()
        return fig

    def plot_custom_feature_box(self, feature, y_title="feature"):
    
        n_clusters = self.clust.chosen_index + 2
        
        colors = sample_colorscale(self.colorscale, np.linspace(0, 1, self.clust.chosen_index+2))

        y_list = []
        for i in range(n_clusters):
            y = feature[np.argwhere(self.clust.labels[self.clust.chosen_index]==i).flatten()]
            y_list.append(y)

        fig = go.Figure()
        for i, y in enumerate(y_list):
            fig.add_trace(go.Box(y=y, marker={"color":colors[i]}, name="Cluster #{}".format(i)))
            
        fig.update_layout(yaxis_title=y_title)
        fig.show()
        return fig

        
    def plot_frames_from_cluster(self,
                                cluster_index, n_frames, 
                                plot_representative=False, 
                                representation="cartoon",
                                extract_selected=False,
                                folder_path="./extracted_frames"):
        
        clust = self.clust
        engen = self.engen
        # Sample frames
        n_clusters = clust.chosen_index + 2
        colors = sample_colorscale("viridis", 
            np.linspace(0, 1, clust.chosen_index+2), colortype="rgb")
        color = colors[cluster_index]
        cluster_frames = np.argwhere(clust.labels[clust.chosen_index]==cluster_index)[:, 0].flatten()
        cluster_size = cluster_frames.shape[0]
        print("Cluster {} contains {} frames/structures".format(cluster_index, cluster_size))
        if n_frames > cluster_size:
            n_frames = cluster_size
        print("Uniformly sampled {} frames from cluster {}:".format(n_frames, cluster_index))
        selected_frames = cluster_frames[np.linspace(0, cluster_size-1, n_frames, dtype=int)]
        print(selected_frames)
        
        # Plot frames with NGLViewer
        frames_mdtraj = []
        nw = nv.NGLWidget()
        for i, frame in enumerate(selected_frames):
            print("Adding frame {}".format(frame))
            traj_frame = mdtraj.load_frame(engen.traj, top=engen.ref, index=frame)
            frames_mdtraj.append(traj_frame)
            nw.add_trajectory(traj_frame,  default_representation=False, name = "Frame {}".format(frame))
            nw[i].add_representation(representation, color=color)
        
        if plot_representative:
            # Plot representative frames in red
            repres_index = clust.chosen_frames[cluster_index]
            print("Adding the representative frame {}".format(repres_index))
            traj_frame = mdtraj.load_frame(engen.traj, top=engen.ref, index=repres_index)
            frames_mdtraj.append(traj_frame)
            nw.add_trajectory(traj_frame,  default_representation=False, name = "Representative frame {}".format(repres_index))
            nw[i+1].add_representation(representation, color="red")
            
        nw.center()
        
        # Extract frames into a separate folder
        if extract_selected:
            # Check if the folder exists
            if not os.path.exists(folder_path):
                # Create the folder
                os.makedirs(folder_path)
            for i, frame in enumerate(selected_frames):
                file_name = os.path.join(folder_path, "cluster_{}_frame_{}.pdb".format(cluster_index, frame))
                print("Extracting {}".format(file_name))
                if os.path.exists(file_name):
                    print("File with the same name exists")
                    print("Deleting the old file")
                    os.remove(file_name)
                
                ret_val = call(["mdconvert -t " + engen.full_ref + " -o " + file_name + " -i " + str(frame) + " " + 
                                engen.full_traj_name], shell=True)
                if not ret_val == 0:
                    raise(Exception("Error saving conformations."))
                    
            # Extract representative too
            repres_index = clust.chosen_frames[cluster_index]
            file_name = os.path.join(folder_path, "cluster_{}_representative_frame_{}.pdb".format(cluster_index, 
                                                                                                repres_index))
            print("Extracting {}".format(file_name))
            if os.path.exists(file_name):
                print("File with the same name exists")
                print("Deleting the old file")
                os.remove(file_name)
            
            ret_val = call(["mdconvert -t " + engen.full_ref + " -o " + file_name + " -i " + str(frame) + " " + 
                            engen.full_traj_name], shell=True)
            if not ret_val == 0:
                raise(Exception("Error saving conformations."))
            
        return nw


        # to avoid loading the full array into memory
    def load_as_memmap(self, pyemma_source):
        pyemma_data_dimension = pyemma_source.dimension()
        pyemma_data_n_frames = pyemma_source.n_frames_total()
        pyemma_iter = pyemma_source.iterator()
        # initialize array
        data_memmap = np.memmap('tmp_data.mymemmap', 
                                dtype='float32', 
                                mode='w+', 
                                shape=(pyemma_data_n_frames,pyemma_data_dimension))
        shape_cnt = 0
        for i, elem in tqdm.tqdm(enumerate(pyemma_iter),total = pyemma_source.n_chunks(pyemma_source.chunksize)):
            data = elem[1]
            data_memmap[shape_cnt : shape_cnt+data.shape[0], :] = data
            shape_cnt += data.shape[0]

        if not shape_cnt == pyemma_data_n_frames:
            print("Warning: output data dimensions don't match the number of frames")
        
        return data_memmap