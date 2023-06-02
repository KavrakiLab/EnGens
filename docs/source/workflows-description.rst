Workflow Description
===================================

.. _workflows-description:

Flow Chart
------------

Here is the visual flow-chart of the pipeline

|EnGens Workflows|
|EnGens Big Logo|

.. |EnGens Workflows| image:: _static/EnGensWF.svg
  :width: 600
  :alt: EnGens Workflows
  
.. |EnGens Big Logo| image:: _static/large-logo.svg
  :width: 600
  :alt: EnGens Big Logo

Workflows detailed description
------------

Workflow 1S: Extracting featurized representations from raw data (static use-case)
------------
Workflow 1S-1: input description
~~~~~~~~~~

The input for the first workflow in the static use-case can be given as any of the following:
List of files in the PDB format 
List of PDB codes (corresponding PDB files are downloaded using pdb-tools package)
Uniprot IDs of a protein (Protein Data Bank is queried for the PDB codes associated the protein with each given 
Uniprot ID using the PDB 1D Coordinate Server; the corresponding PDB files are fetched using pdb-tools package) 


All of these inputs effectively result in a list of files in the PDB format, which contain Cartesian coordinates of 
all atoms in the protein structure. Before proceeding with the downstream workflows these input structures are preprocessed.

Workflow 1S-2: preprocessing input
~~~~~~~~~~

Structures are first cleaned with the `PDBFixer`_ tool. By default the preprocessing includes: 
* adding missing residues and atoms; 
* removing ligands, hydrogens, water and other heteroatoms. 
Users can edit these default settings.

Next, it is necessary to renumber the residues listed in these PDB files, so that each file involves the same numbering. 
Residues are renumbered using `PDBrenum`_ tool according to their UniProt sequence.
PDB files commonly contain multiple chains and this requires additional care during preprocessing steps. Often, these multiple chains correspond 
to different domains or subunits of the protein (e.g., catalytic and regularly subunits). In order to process these domains we fetch the 
metadata related to the PDB codes from the Protein DataBank using the `RCSB-API`_. This metadata provides information such as Uniprot IDs and protein names of 
each included chain. Users can look at these values and decide which chains to process and if they want to discard certain chains. 
The metadata might look something like this:

|PDB metadata|

Downstream analysis is performed on all chains (asymetric ids) selected by the user. 
Note that experimentally-resolved structures of the same protein may vary in the number of resolved residues. 
Some structures may contain additional 
chains and artifacts that helped stabilize the protein for crystallization. 
thers may lack regions that were not resolved. This is undesirable for 
downstream structural analysis. We want to analyze only the relevant parts of the structure and we want every structure to have the same number of residues. 
For these reasons, we first have to “trim” the structures to the same number of residues. To this end, we perform multiple structure alignment of the input 
files using `mTM-align`_ for multiple structure alignment. As part of the structure alignment, mTM-align provides a sequence mapping where the amino acids of the sequences are matched 
to each other based on the structure alignment - we identify continuous regions of these amino acids (where no structure has a gap). 
These continuous regions form what we call the maximum common substructure (MCS). We extract the residues corresponding to the MCS from each structure. As a result, 
the input structures are “trimmed” to an MCS with a constant number of residues. Users can visually inspect this sequence 
mapping produced by the structure alignment and the extracted substructures to make sure that they are valid. Example of the alignment output and 
the MCS structure are shown bellow

|mTM alignment|
|MCS structure|

Note that the manual 
inspection of the MCS regions from the alignment is an extremely important step. If the insufficient part of the protein structure is extracted with MCS, 
the analysis might be performed without an important structural component.
Additionally, after extracting the MCS, users might be interested in analyzing only a specific part of a protein (e.g. only the binding pocket). 
We provide an option for selecting a region by identifying the residues of interest (using the selection algebra of the mdtraj package). 
With this, the preprocessing of input structures is finalized and all the downstream analyses are performed on this preprocessed region of the protein.

.. _PDBFixer: https://github.com/openmm/pdbfixer
.. _PDBrenum: http://dunbrack.fccc.edu/PDBrenum/
.. _RCSB-API: https://data.rcsb.org/#gql-api
.. _mTM-align: https://yanglab.nankai.edu.cn/mTM-align/

.. |PDB metadata| image:: _static/pdb_metadata.jpg
  :width: 1000
  :alt: PDB metadata
  
.. |mTM alignment| image:: _static/mTM-align.jpg
  :width: 800
  :alt: mTM alignment
  
.. |MCS structure| image:: _static/MCS_structure.JPG
  :width: 400
  :alt: MCS structure
  
Workflow 1S-3: featurization
~~~~~~~~~~

Each preprocessed input structure :math:`\s` is still represented by the x, y, z Cartesian coordinates of all its  N  atoms - i.e.,  :math:`s \\in R^3`. 
In other words, structure s represents a point in configuration space R3N. This representation is redundant - 
the motion of atoms is constrained by bonds and the underlying configuration space of the structure is actually a manifold of :math:`R^3`.  
Additionally, this representation is not invariant to translation and rotation. 
Various mappings can be applied to embed each input structure into a different, more desirable space. 
Such mapping is often referred to as featurization.

For example, this mapping can be applied to extract  and  torsion angles of the protein backbone. 
In this case  has the non-Euclidian topology of a r-dimensional torus :math:`T^r` (where r is the number of torsion angles). 
Note that standard algorithms for dimensionality reduction and clustering perform statistical analysis on Euclidean data, 
as Euclidean distance is used when calculating the mean and covariance for example. Approaches for dealing with non-Euclidean 
manifolds have been proposed, but are out of the scope of our current work. Thus,  must be Euclidean to comply with the downstream algorithms. 
A simple work-around recommended for dealing with angular values is to map angles to their sine and cosine values [sin(), cos()]. In this case,   
would effectively be R2r.  
Mapping  can be tailored to a specific protein and motivated by a biological insight. 
For example, it is known that the relative position of the nSH2 domain with respect to other domains of 
PI3K distinguishes active and inactive states of this enzyme. In this case an interesting featurization would contain 
pairwise distances between PI3K domains:  can featurize the structure by calculating pairwise distances between the centers of mass 
of these domains. In that case  would correspond to Rd(d-1)2 (where d is the number of domains). 
In both examples, the resulting featurized representations produced by  will be invariant 
to rotations and translations of the structure in the original space. 
We use the implementation of featurizations provided by the `PyEmma`_ package which includes multiple featurization options summarized in Table S2. 
Note that the choice of featurization can be challenging and highly depends on the system in question. 
For example, when studying conformational changes in the binding site of a protein, 
the user could consider using residue-level distance-based features within the binding pocket. On the other hand, 
if it is necessary to pay attention to the distance between substrates and catalytic residues, 
the user can choose distance-based features on an atomistic level. When studying big conformational changes of the whole protein, 
larger portions (domains) of the protein should be considered. In that case, it is common to use flexible torsions of the whole backbone, 
and distance-based features for C-alpha atoms or groups of residues.  When studying the interaction between two domains (or a loop and a domain), 
one good choice is to use pairwise distances between the centers of mass (COM) of the domains. 
User input for the featurization is very important and it is advisable to test out a couple of options and inspect the results for consistency. 

.. _PyEmma: http://emma-project.org/latest/

Workflow 1S-4: output description
~~~~~~~~~~

After applying the featurization  to each structure in the dataset, 
we combine all structures into a matrix :math:`D_{nm}`, where  n is the number of input structures and m is the dimensionality of the featurization. 
Note that m is usually very high. For example, if the center of mass of the backbone C-alpha atoms is used for featurization of a protein with 300 residues, 
the resulting dimensionality m will be 900. This motivates the second workflow, dimensionality reduction.


Workflow 1D: Extracting featurized representations from raw data (dynamic use-case)
------------
Workflow 1D-1: input description
~~~~~~~~~~
In the dynamic use-case the input is the result of a molecular dynamics simulation, which includes the following:
a trajectory (in any format accepted by `PyEmma`_: .xtc, .dcd, .hdf5, .trr)
a topology (in PDB format)

Workflow 1D-2: preprocessing input
~~~~~~~~~~

Frames in the trajectory should first be aligned using mdtraj package to a reference frame. 
This is especially important when features dependent on translation and rotation 
(e.g., any featurization expressed in atomic Cartesian coordinates) are used, since small movements of the whole protein can add a lot of noise.
The next step involves selecting a region of interest (such as a binding site). 
Note that, contrary to the static use-case, here it is not necessary to trim the structures to an MCS - 
all the frames have exactly the same topology and the same number of residues.

.. _mdtraj: https://www.mdtraj.org/1.9.8.dev0/index.html

Workflow 1D-3: featurization
~~~~~~~~~~

By default, each frame is represented by the Cartesian coordinates of all atoms. 
Similarly to the static use-case, featurization (see Workflow 1S-3) can be viewed as a transformation :R3N  
applied to each frame. Frames can be better represented by any of the featurizations provided by PyEmma. 
Few related works have tried to address the issue of featurization. 
These works aimed to devise general advice on choosing a proper featurization regardless of the system in question.  
Note that most of these works address the question of featurization in context of building Markov State Models from MD trajectories. 
We assume that many of their findings still apply in the context of generating representative conformational ensembles aith EnGens.
Husic et al. (2016) investigated the performance of using whole-backbone features including: α-angles, α-carbon contact distances, pairwise α-carbon RMSD. 
They did not report any consistent trends in featurization. Scherer et al. (2019) investigated the use of aligned cartesian coordinates, distance-based features, 
contact features, solvent accessible surface area, flexible torsions, and combined distance transformations and flexible torsions. 
They found that using torsions alone leads to inferior models but that combining flexible torsions with contact information improves models.
Studying conformational changes in biomolecules is a very complex task. 
It is advisable to try out and compare different featurization schemes. User input is still essential to understand which featurization performs 
best for the phenomena of interest. 


Workflow 1D-4: VAMP featurization tuning
~~~~~~~~~~

Choosing the featurization is very important in both use-cases, as it can greatly impact the downstream analysis. 
In general, an optimal featurization will allow for identifying biologically-relevant and interesting states. One 
way to formalize this is to say that biologically relevant events are rare and that the processes that govern them 
are slow. A good featurization would allow for the identification of such slow processes. 
Recent advances in Markov modeling of the system dynamics from MD data provide insights into resolving slow processes 
that guide the dynamics. Many approaches have been developed in an effort to model the dynamics from MD data, including 
time invariant component analysis (TICA), markov state modeling (MSM), variational approach to conformational dynamics (VAC), 
variational approach for learning Markov processes (VAMP). These methods are able to identify dominant (slow) dynamics by 
approximating the transfer operator guiding the dynamics. The VAMP approach approximates the Koopman transfer operator of the
underlying Markovian dynamics by deriving its top singular components from the time series data. This approach is valuable 
because it can deal with both reversible and irreversible dynamics and because it provides VAMP-scores. The VAMP-E score 
quantifies the error in the approximated Koopman operator and can be used to tune the VAMP hyper-parameters. One of these 
hyper-parameters is the featurization. Thus, the featurization that yields the maximum VAMP-E score allows for the most 
accurate resolution of the dynamics and can be viewed as optimal.
Here, we use PyEmma’s implementation of VAMP scores to allow users to evaluate their featurization by calculating 
the VAMP-E score. Users can choose from a range of different featurizations by picking the one with the highest VAMP score.

  Hao Wu and Frank Noé. Variational approach for learning markov processes from time series data. Journal of Nonlinear Science, 30(1):23–66, 2020.
  
  Frank Noé and Cecilia Clementi. Kinetic distance and kinetic maps from molecular dynamics simulation. Journal of chemical theory and computation, 11(10):5002–5011, 2015.
  
Workflow 1D-5: output description
~~~~~~~~~~
After applying the selected featurization to each frame in the trajectory - frames are combined into a matrix 
:math:`D_{nm}`, where  n is the number of frames and m is the dimensionality of the featurization. 
This matrix is similar to the matrix resulting from Workflow 1S and almost all downstream analyses will 
be the same for the static use-case (S) and the dynamic use-case (D). The only exceptions are dimensionality 
reduction techniques in steps 2D-2c and 2D-2d (i.e., TICA and SRV). 
These methods model the dynamics of the data and can thus be applied only to the dynamic use-case.


Workflow 2(S&D):  Projecting the featurized representation into an embedding in low-dimensional space (static and dynamic use-case)
------------

Workflow 2-1: input description
~~~~~~~~~~

The input to Workflow 2 consists of the featurized dataset expressed as the matrix :math:`D_{nm}`. 
For the static use-case n corresponds to the number of collected structures, for the dynamic use-case n is the number of frames in the trajectory. 
m corresponds to the dimensionality of the chosen featurization. Since m is usually high as the featurizations have high dimensionality, 
the important step executed in Workflow 2 is reducing this dimensionality by projecting the data into a lower dimensional 
space by preserving important properties of the data expressed as the variance or the kinetic variance. 
Different methods offered by EnGens are presented in following subsections including PCA (WF2-2a), TICA (WF2-2b), UMAP (WF2D-2c) and SRV (WF2D-2d). 
Note that PCA and TICA are linear methods that find a linear projection into a lower-dimensional space, while UMAP and SRV perform a non-linear projection. 
In addition, TICA and SRV use the dynamic properties of the underlying MD datasets and should not be applied to the static use-case.

Workflow 2-2a: PCA
~~~~~~~~~~
Principal component analysis (PCA) is a widely used and powerful statistical technique for reducing dimensionality. 
The key idea of PCA is in finding the “principal components” (PCs - vectors in high-dimensional spaces) onto which the data 
can be projected while preserving as much variance as possible. 
The first step is to normalize the data by subtracting the mean and the standard deviation of the dataset from each point. 
This is important to ensure that individual feature variances are comparable. Next, the eigenvectors of the data metrics are calculated using  
singular value decomposition methods. These eigenvectors correspond to PCs. Next, the data is projected onto the PCs. The k PCs with the highest 
eigenvalues are selected, effectively reducing the dimensionality of the input data to k dimensions. Usually, k is chosen such that a satisfactory 
amount of variance is retained (e.g., 90%).
EnGens uses scikit-learn implementation of the PCA method to obtain the final matrix :math:`D_{nk}` where n is the number of structures (or frames) 
in the dataset and k is the number of chosen PCs. Note that PCA is a linear method: it assumes that the relationship between features can be represented 
by a linear combination of the original features. If a nonlinear relationship between features is suspected, a more suitable method for 
dimensionality reduction might be UMAP.

  Tipping, M. E., and Bishop, C. M. (1999). “Probabilistic principal component analysis”. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 61(3), 611-622

Workflow 2-2b: UMAP
~~~~~~~~~~
The Uniform Manifold Approximation and Projection (UMAP) algorithm is a non-linear dimensionality reduction technique that can preserve 
the global structure of the data. UMAP aims to approximate the underlying manifold of the data by constructing a fuzzy topological 
representation of these data and projecting them onto a lower-dimensional space.
The first step in UMAP is the construction of the weighted fuzzy graph. This step is controlled by parameters: n_neighbors and min_dist. 
These parameters indicate the number of nearest neighbors that each point is connected to and the minimum distance between neighbors. 
Next, the spectral embedding of the constructed graph is calculated and optimized with respect to the fuzzy set cross-entropy 
(or its differentiable approximation) using stochastic gradient descent. The terms of the fuzzy set cross-entropy objective function 
can be viewed as attractive and repulsive forces. Optimizing this objective will balance the “pull” and “push” of the data and the 
low-dimensional representation should settle into a state that relatively accurately represents the overall topology of the original data. 
The n_components parameter indicates the dimensionality of the embedding and it has to be selected prior to fitting UMAP. 
EnGens uses an open-source UMAP implementation to produce the final matrix :math:`D_{nk}` where n is the number of structures (or frames) 
in the dataset and k is equal to n_components (the number of chosen components). 
Note that UMAP assumes that the data are uniformly distributed on the underlying manifold. 
Consequently, it does not preserve the relative local density of the data well. 
This fact is highlighted when the two-dimensional free energy plots of UMAP projections are generated, 
as they rely on estimating the density of clusters.  
Thus, when using UMAP for dimensionality reduction, all downstream analyses of density estimation should be performed with care.

  McInnes et al., (2018). UMAP: Uniform Manifold Approximation and Projection. Journal of Open Source Software, 3(29), 861

Workflow 2D-2c: TICA
~~~~~~~~~~
Time-invariant component analysis (TICA) is a dimensionality reduction technique commonly used in the analysis of time series data 
(such as MD trajectory data). Similarly to PCA, TICA finds the components (vectors) onto which all datapoints should be projected. 
However, the objective of TICA is not to find the components preserving most of the variance. Instead, TICA aims to identify slow 
dynamics that govern the behavior of the system. To this end, TICA looks for a coordinate system by identifying the linear combinations 
of the original coordinates that change most slowly over time.
The optimization problem formalized by TICA is to find coordinates of maximal autocorrelation at the given lag time. This can be achieved 
by the following few steps. First, the  mean-free covariance (C0) and mean-free  time-lagged covariance (C) matrices are calculated from the 
data. Next, the generalized eigenvalue problem is solved Cri =C0iri , where riare the eigenvectors corresponding to time independent components
- TICs,  while i are the respective eigenvalues that indicate the relaxation time-scales of the TICs. Finally, the k slowest resolved processes 
(TICs) are picked to project the data. Choosing k can be done  by setting a threshold on the kinetic variance that should be retained in the data. 
EnGens uses PyEmma’s implementation of TICA to produce the final matrix Dnk where n is the number of frames in the trajectory and k is the number 
of chosen TICs. Note that TICA is a linear method and assumes a linear relationship between features. For resolving non-linear slow processes the 
SRV method is recommended.

  Perez-Hernandez, G. et al. (2013). Identification of slow molecularorder parameters for Markov model construction. The Journal of Chemical Physics, 139(1), 015102
  
  Schwantes, C. R. and Pande, V. S. (2015). Modeling Molecular Kinetics with tICA and the Kernel Trick. Journal of Chemical Theory and Computation, 11(2), 600–608
  
Workflow 2D-2d: SRV
~~~~~~~~~~
The state-free reversible VAMPnet (SRV) is a non-linear dimensionality reduction technique. 
Similarly to TICA, SRV aims to identify the slow processes (modes) that drive the dynamics of the system. 
In contrast to TICA, which solves a generalized eigenvalue problem, SRV resolves the processes by learning the leading 
eigenfunctions of a dynamical system transfer operator from trajectory data. Using a deep learning architecture allows SRV to 
identify non-linear processes.
In SRV, a Siamese neural network is constructed with two subnets that share the same architecture and weights. A pair of 
datapoints is fed into the network. The pairs are selected to be rows of the featurized matrix Dnmseparated by the lag time (dt, dt+:). 
The network learns the mappings fi, where i=(1,..,k) and k is the number of modes we want to learn. This parameter is set by the user; 
it defines the architecture as the number of nodes in the last hidden layer of the network will be k. The network produces the outputs (fi(dt), fi(dt+:)). 
The loss function is calculated using these outputs and it corresponds to the negative VAMP-2 score. 
As the network maximizes the VAMP-2 score, it learns non-linear feature embeddings that allow for the best approximation of the transfer operator. 
In effect, this results in learning a good non-linear low-dimensional embedding of the original features.
EnGens uses the open-source implementation of SRV  to produce the final matrix Dnk where n is the number of frames in the trajectory and k is 
the selected number of modes. When SRV is used, the embedding dimensionality k and the lag time  have to be carefully selected by the user. 

  Chen, W. et al. (2019). Nonlinear discovery of slow molecular modes using state-free reversible VAMPnets. The Journal of Chemical Physics, 150(21), 214114

Workflow 2-3: output description
~~~~~~~~~~

After projecting the featurized representation of each datapoint to a lower-dimensional embedding, 
the resulting data matrix has the form :math:`D_{nk}`, where  n is the number of input structures or the number of frames in the trajectory, 
and k is the dimensionality of the featurization. 

Workflow 3(S&D): Clustering embeddings and extracting the ensemble (static and dynamic use-case)
------------

Workflow 3-1: input description
~~~~~~~~~~
For both the dynamic and static use-case, the input to Workflow 3 consists of the embedded dataset expressed as the matrix :math:`D_{nk}`. 
For the static use-case n corresponds to the number of collected structures, for the dynamic use-case n is the number of frames in the trajectory. 
k corresponds to the dimensionality of the embedding. After reducing the dimensionality of the data in Workflow 2 it is possible to visualize the dataset 
in a 2D or 3D space defined by the top identified component, which gives better insight into the dataset. 
The additional step of clustering performed in Workflow 3 allows for the data to be stratified into groups of datapoints (structures) 
that are close together in this embedded space and are assumed to be structurally similar. 
To this end, users can choose one of these widely used methods: K-means, hierarchical or GMM clustering.

Workflow 3-2a:  K-means clustering
~~~~~~~~~~
The K-means algorithm aims to separate n datapoints into c clusters by assigning each datapoint to one of the c centroids. 
First, the positions of the centroids are initialized (randomly or via some heuristic). Next, the algorithm iterates between: 
(1) generating clusters by assigning points to the closest centroid, and (2) repositioning the centroids to the mean of the clusters. 
The centroid assignment is evaluated using inertia (or sum of distances within a cluster).  
The algorithm iterates until the movement of the centroids between subsequent steps is under some threshold (i.e., until convergence). 
The number of centroids c (i.e., the number of produced clusters) has to be predefined. 
This value is often difficult to predict and it is recommended users evaluate different values. 
After running K-means with a range of c values, the best c can be selected by plotting the corresponding inertia and 
using the elbow method (finding the c for which inertia drops rapidly). Using the elbow method, c will be selected such that 
clusters are compact around the centroids. Alternatively, users can inspect the silhouette plots for different c values and 
select the value with the maximum average silhouette score. The silhouette score is assigned to each datapoint.
The silhouette score is in the range (-1, +1) where -1 indicates overlapped clusters and +1 indicates tight, well separated clusters.  
Maximizing the silhouette score will result in picking the clustering with clear separation between datapoints and little overlap between clusters.
EnGens uses `scikit-learn`_ implementation of the K-means algorithm to generate the cluster assignment of all datapoints. 
Note that K-means uses inertia to evaluate clusters, which relies on convex and isotropic assumptions. 
This results in spherical clusters, which might not always be a good fit for the data. 
Alternatively, if the underlying data distribution is Gaussian, GMM clustering can be used to produce clusters of different shapes, sizes, and orientations.

  Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A K-Means Clustering Algorithm. Journal of the Royal Statistical Society. Series C (Applied Statistics), 28(1), 100–108.

.. _scikit-learn: https://scikit-learn.org/stable/

Workflow 3-2b: GMM clustering
~~~~~~~~~~
Gaussian Mixture Models are probabilistic generative models that can represent a dataset as a mixture of Gaussian distributions with unknown parameters. 
The number of distributions (components) is predefined (similarly to the number of centroids in K-means). 
The expectation-maximization (EM) algorithm is employed to fit the parameters of the Gaussian distributions to the data. 
EM consists of consecutive steps of (E) calculating the log-likelihood of the data given the current set of parameters and 
(M) finding new parameters that maximize the given log-likelihood.  After modeling the data as a mixture of Gaussian distributions, 
each datapoint ends up associated  with a probability for each of the GMM components. Eventually, 
datapoints are assigned to the GMM component (i.e., the cluster) with the highest associated probability. 
The number of GMM components is predefined. Users can optimize this number by running GMM on a range of values and picking the one with the 
lowest associated Bayesian information criterion (BIC). BIC is defined as k ln(N) - 2ln(L), where k is the number of GMM components, 
N is the number of datapoints and L is the likelihood of the data under the inferred GMM model. BIC penalizes model complexity with the first 
term (preventing overfitting) and rewards model fitting with the second term. Analysis of BIC allows picking k using the elbow method 
(finding the number of components for which BIC drops rapidly). Alternatively, clustering can be assessed with the silhouette method in a 
similar manner as for the K-means method.
EnGens uses `scikit-learn`_ implementation of the GMM algorithm to produce the cluster assignment of the datapoints. 
Note that both K-means and GMM require careful consideration of the number of resulting clusters, and that (a range of values for) 
this parameter has to be defined before running the algorithm. On the other hand, hierarchical clustering uses a different approach 
which allows defining the number of clusters after computing a dendrogram of the data. 

  Lindsay, B. et al. (1989). Mixture Models: Inference and Applications to Clustering. In Journal of the American Statistical Association, volume 84, page 337.

Workflow 3-2c: hierarchical clustering
~~~~~~~~~~
Hierarchical clustering forms groups of datapoints by successively merging or dividing the dataset. 
In EnGens, we employ the agglomerative hierarchical clustering, which is a bottom-up approach which iteratively merges datapoints into groups. 
First, each datapoint defines its own cluster. Next, clusters are iteratively merged based on a given criterion. For example, when the average 
linkage criterion is used, a pair of clusters is merged if the average distances between all their points are minimum. 
Merging steps are performed until only one cluster remains. The resulting dendrogram can be inspected to choose the number of clusters 
that best describes the data. 
EnGens uses `scikit-learn`_ implementation of the agglomerative hierarchical clustering. Note that for a large number of datapoints, 
agglomerative clustering can be computationally intensive, which might be of particular concern for MD data with a high number of frames. 
To combat this issue, users can provide a subsample of the trajectory as input (e.g., by loading every 10th frame) by setting the stride parameter. 

  Murtagh, F. and Contreras, P. (2012). Algorithms for hierarchical clustering: an overview. WIREs Data Mining and Knowledge Discovery, 2(1), 86–97. 

Workflow 3-3: representative selection
~~~~~~~~~~

After clustering, the input matrix :math:`D_{nk}` , each cluster contains a set of datapoints (structures) that are assumed to be structurally similar. 
In Workflow 3-3,  a single datapoint, called the cluster representative, is selected from each cluster to represent it within the final representative ensemble. 
If users do not want to pick representatives from all clusters (e.g., to ignore very small clusters), a weight threshold can be used to discard clusters.
There are two ways to pick a representative from each cluster. First, the representative can be the point closest to the center of the cluster: 
in K-means, cluster centroids are the centers; in GMM, the means of the distributions are the centers.
Alternatively, the representative can be a hub of the cluster, which is defined as the point with the most neighbors in the cluster.
Neighbors are defined as datapoints that are closer to the hub than the mean of pairwise distances within the cluster.
   
Workflow 3-3: output description
~~~~~~~~~~
The result of this workflow is (1) the grouping of the datapoints into clusters of structures with similar properties and (2) the 
final representative ensemble composed of the extracted cluster representatives. 

Workflow 4(S&D):  Visualizing the data and analyzing the ensemble
------------
EnGens provides a variety of plots to summarize the resulting data clustering and generated representative ensemble. 
Plots are implemented using `plotly`_, `NGLViewer`_, `py3Dmol`_ is used to visualize the 3D structures of the representatives. 

.. _plotly: https://plotly.com/
.. _NGLViewer: https://nglviewer.org/
.. _py3Dmol: https://github.com/avirshup/py3dmol

Embedding plot
~~~~~~~~~~
This plot shows the projected embedding from Workflow 2 in two dimensions. 
The X-axis represents the first principal component (or TIC/SRV/UMAP component) while the Y-axis represents the second principal component 
(or TIC/SRV/UMAP component). Each point corresponds to one structure (or frame) from the dataset. Points are colored based on their cluster membership.

Cluster membership plot
~~~~~~~~~~
This plot shows the details of cluster membership across the whole dataset. 
Y-axis represents cluster membership. 
In the static case, the X-axis lists the structure IDs in the order specified in the original dataset. 
In the dynamic case, the X-axis lists the MD frames according to the simulation time.

Multiple structures visualization
~~~~~~~~~~
We use NGLViewer to show the extracted representative ensemble within the notebook. 
It is possible to select a cluster of interest and display the 3D cartoon representation of the cluster representative. 
In addition, other uniformly sampled structures from the same cluster can be visualized. 

Ensemble summary box plots
~~~~~~~~~~
This plot summarizes a property of interest  (e.g., RMSD to a reference structure) across different clusters with a box plot. 
The X-axis lists the cluster IDs. The Y-axis corresponds to the selected property. 

Ensemble summary scatter plots
~~~~~~~~~~
This plot summarizes a property of interest  (e.g., RMSD to a reference structure) across different clusters with a scatterplot. 
The X-axis lists the structure IDs (i.e, the number of an MD frame in the dynamic data, or a structure name in the static data). 
The Y-axis corresponds to the selected property. Points are colored based on their cluster membership. 
Selected representative structures are highlighted in red.
