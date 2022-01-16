from logging import FileHandler
import sys
import time
import mysql.connector
from subprocess import Popen
from config import *
import os
sys.path.append(os.path.abspath('../engens_code/'))
from engens.core.EnGens import EnGen
import engens.core.FeatureSelector as fs
from engens.core.DimReduction import *
from engens.core.ClustEn import *
import pickle as pk
import mdshare
import nglview as nv
from outputWriter import *

cnx = mysql.connector.connect(user = 'root', password = 'admin12345678', host = '127.0.0.1', database = 'kavraki')
recordId = sys.argv[1]
cur = cnx.cursor()
sqlQuery = "SELECT * FROM uploads WHERE resultId = %s;"
cur.execute(sqlQuery,[recordId])

referencePDBfilePath = ""
referenceXTCfilePath = ""
for (id, pdb_file_name, upload_time, xtc_file_name, access, pass_salt, password, username, resultId, title, email) in cur:
    referencePDBfilePath = upload_folder + pdb_file_name
    referenceXTCfilePath = upload_folder + xtc_file_name
    print("Email: ", email, "Project title:", title, "Access", access)

st = outputWriter(resultId)
st.title("Ensemble Generation Notebooks (EnGeNs)")
st.write("Email: " + email) 
st.write("Project title: " + title) 

cur.close()
cnx.close()

print("PDB file: " + referencePDBfilePath)
print("XTC file: " + referenceXTCfilePath)

pdb = mdshare.fetch('pentapeptide-impl-solv.pdb', working_directory='.')
files = mdshare.fetch('pentapeptide-00-500ns-impl-solv.xtc', working_directory='.')

st.wf1intro()
st.wf1step1()

#location of trajectory and topology files
top_loc = referencePDBfilePath
traj_loc = referenceXTCfilePath
engen = None

if engen_object_config == 1:
    # Example 1 - load full trajectory from the files and visualize it
    #instrantiate the engen object
    engen = EnGen(traj_loc, top_loc)

if engen_object_config == 2:
    # Example 2 - load a subset of trajectory using a list of atom indices
    #for example - first N atoms?
    N = first_N_atoms
    select_list = [i for i in range(N)]
    #instrantiate the engen object
    engen = EnGen(traj_loc, top_loc, select_list)

if engen_object_config == 3:
    # Example 3 - load a subset of trajectory using the atom selection string

    #build the selection string
    binding_site_selstr = ""
    for res_ind in binding_site_residues:
        if binding_site_selstr=="": 
            binding_site_selstr = "resid "+str(res_ind)
        else:
            binding_site_selstr += " or resid "+str(res_ind)

    print("Resulting binding site selection string:")
    print(binding_site_selstr)
    st.write("Resulting binding site selection string:" + binding_site_selstr)

    #instrantiate the engen object
    engen = EnGen(traj_loc, top_loc, binding_site_selstr)

#visualize the trajectory (optional - if trajectory too large, skip this step)
nglwidget = engen.show_animated_traj()
#st.writeNGL(nglwidget, "ngltrajectory")

st.wf1step2()

# initialize default features 
engen.init_featurizers_default()

# print the desctiption of the default features
description = engen.describe_featurizers()
print(description)
#st.write(description)

# make a selection of functions to be applied from pyemma featurizers
# choose from http://www.emma-project.org/latest/api/generated/pyemma.coordinates.featurizer.html
# with respective parameters

#residue mindist
feat1 = {
    "add_all": {}
}
#add the respective features to the engen structure
engen.add_featurizer(feat1)

#C-alpha distances
feat2 = {
    "add_all": {},
    "add_distances_ca": {"periodic":True, "excluded_neighbors":3}
}
#add the respective features to the engen structure
engen.add_featurizer(feat2)

#center of mass and torsion angles
feat3 = {
    "add_residue_COM": {"residue_indices": [1,2,3]},
    "add_backbone_torsions": {"cossin":True, "periodic":False}
}
#add the respective features to the engen structure
engen.add_featurizer(feat3)

# print the desctiption of the selected features features
description = engen.describe_featurizers()
print(description)
st.write(description)

st.wf1step3()

#initialize some features - default of custom
engen.init_featurizers_default()

#provide a range of lags and dimensions to train on
lags = [100, 200, 300, 500, 1000, 1500]
dims = [i + 1 for i in range(1,100, 20)]

#initialize VAMP2 featurizer and run it
featsel = fs.VAMP2FeatureSelection(lags, dims, engen)
#--featsel.run_vamp()

#plot VAMP2 results
#--featsel.plot_results()

#plot how number of dimensions increases VAMP2 score for different lag times for a specific featurization

feat_index = 2
#fig = featsel.plot_dimensions(feat_index)
#st.write("Plot how number of dimensions increases VAMP2 score for different lag times for a specific featurization")
#st.writeFig(fig, "featsel_index_2_")


st.wf1step4()

# Option 1 - select using the VAMP2 selector from above 
#--featsel.select_feature()

# Option 2 - select using your own analysis, just set the number of the feature from the list
feat_num = 1

# initialize selector
featsel = fs.UserFeatureSelection(feat_num, engen)
#select the feature
feature = featsel.select_feature()
st.write(feature)

st.wf1step5()

# ADDED BY ALEX
engen.apply_featurizations()

with open("./static/outputfiles/" + resultId + "_wf1_resulting_EnGen.pickle", "wb") as file:
    pk.dump(engen, file, -1)

print("WF1 finished. Generated file:" + "static/outputfiles/" + resultId + "_wf1_resulting_EnGen.pickle")

st.wf2intro()

#with open("static/outputfiles/" + resultId + "_wf1_resulting_EnGen.pickle", "rb") as file:
#    engen = pk.load(file)

st.wf2step1()

traj = engen.traj
ref = engen.ref
print("Using the trajectory {} and reference pdb file {}".format(traj, ref))

topology = engen.mdtrajref
print("The topology is:")
print(topology)
st.write(topology)

print(engen.chosen_feat_index)

frame_num = engen.data[engen.chosen_feat_index][1].shape[0]
print("Number of frames is {}".format(frame_num))
st.write("Number of frames is {}".format(frame_num))

feat_dims = engen.data[engen.chosen_feat_index][1].shape[1]
print("The dimensionality of your featurization is {}".format(feat_dims))
st.write("The dimensionality of your featurization is {}".format(feat_dims))

feat = engen.featurizers[engen.chosen_feat_index]
print("You chose to featurize with")
print(feat.describe())
st.write(feat.describe())

st.wf2step2()

#reduction_type = "TICA" #can be "TICA"|"HDE"|"PCA"
reducer = dimreds[reduction_type](engen)

st.wf2step3()

#analyse the above plot and choose a lag that is closest to the saturation of timescale
reducer.plot_lag_analysis()
#chosen_lag = 10
saveloc = "./static/outputfiles/plots/" + resultId + "_1_plot_lag_analysis.png"
st.writePicture(saveloc)
reducer.plot_lag_analysis(chosen_lag=chosen_lag, save_loc=saveloc)
reducer.choose_lag(lag=chosen_lag)

st.wf2step4()

saveloc = "./static/outputfiles/plots/" + resultId + "_2_plot_2d.png"
st.writePicture(saveloc)
reducer.plot_2d(save_loc=saveloc)

saveloc = "./static/outputfiles/plots/" + resultId + "_3_plot_variance.png"
st.writePicture(saveloc)
if (reduction_type == "TICA"):
    # var_thr = 70
    reducer.plot_variance(var_thr=var_thr, save_loc=saveloc)
    #now you can pick the proposed number of components
    print(reducer.transformed_data.shape)
    #comp_num = 16
    reducer.choose_n(comp_num)
    print(reducer.transformed_data.shape)

st.wf2step5()

# to apply the dimensionality reduction to the engen structure run:
reducer.apply()

# to save the results from this analysis for workflow3 run:
with open("./static/outputfiles/" + resultId + "_wf2_resulting_EnGen.pickle", "wb") as file:
    pk.dump(engen, file, -1)
with open("./static/outputfiles/" + resultId + "_wf2_resulting_Reducer.pickle", "wb") as file:
    pk.dump(reducer, file, -1)

st.wf3intro()

#engen = None
#with open("static/outputfiles/" + resultId + "_wf2_resulting_EnGen.pickle", "rb") as file:
#    engen = pk.load(file)

st.wf3step1()

traj = engen.traj
ref = engen.ref
print("Using the trajectory {} and reference pdb file {}".format(traj, ref))

topology = engen.mdtrajref
print("The topology is:")
print(topology)

frame_num = engen.data[engen.chosen_feat_index][1].shape[0]
print("Number of frames is {}".format(frame_num))

feat_dims = engen.data[engen.chosen_feat_index][1].shape[1]
print("The dimensionality of your featurization is {}".format(feat_dims))

feat = engen.featurizers[engen.chosen_feat_index]
print("You chose to featurize with")
print(feat.describe())

dimred_data = engen.dimred_data
print("After dimensionality reduction the dimension of your features is {}".format(dimred_data.shape[1]))

st.wf3step2()

# Example 1 - choose Kmeans

# clustering = "KM" # options: KM | GMM
cluster_method = clusterings[clustering](engen)

cluster_method1 = clusterings["KM"](engen)
cluster_method2 = clusterings["GMM"](engen)

st.wf3step3()

# Example 1 - iterate through Kmeans number of clusters

# clusterRange = range(2, 10)

params = [{"n_clusters":i} for i in clusterRange]
outputText = cluster_method1.cluster_multiple_params(params)
st.write(outputText)

saveloc = "./static/outputfiles/plots/" + resultId + "_4_elbow_KM.png"
st.writePicture(saveloc)
# analyze these parameters with the elbow method
params, outputText = cluster_method1.analyze_elbow_method(filename=saveloc)
st.write(outputText)

saveloc = "./static/outputfiles/plots/" + resultId + "_5_silhouette_KM.png"
st.writePicture(saveloc)
# analyze these parameters with the silhouette method
outputText = cluster_method1.analyze_silhouette()
st.write(outputText)

# pick the parameter index
# param_index = 2

cluster_method1.choose_param(param_index_KM)

# Example 2- iterate through GMM number of components

params = [{"n_components":i} for i in clusterRange]
outputText = cluster_method2.cluster_multiple_params(params)
st.write(outputText)

saveloc = "./static/outputfiles/plots/" + resultId + "_6_elbow_GMM.png"
st.writePicture(saveloc)
# analyze these parameters with the elbow method
params, outputText = cluster_method2.analyze_elbow_method(filename=saveloc)
st.write(outputText)

saveloc = "./static/outputfiles/plots/" + resultId + "_7_silhouette_GMM.png"
st.writePicture(saveloc)
# analyze these parameters with the silhouette method
outputText = cluster_method2.analyze_silhouette()
st.write(outputText)

# pick the parameter index
# param_index = 0

cluster_method2.choose_param(param_index_GMM)

st.wf3step4()

saveloc = "./static/outputfiles/plots/" + resultId + "_8_cluster_weight_KM.png"
st.writePicture(saveloc)
cluster_method1.plot_cluster_weight(save_loc = saveloc)
saveloc = "./static/outputfiles/plots/" + resultId + "_9_choose_cluster_weight_KM.png"
st.writePicture(saveloc)
outputText = cluster_method1.choose_clusters(thr_KM) # thr_KM = 0.2
st.write(outputText)

saveloc = "./static/outputfiles/plots/" + resultId + "_10_cluster_weight_GMM.png"
st.writePicture(saveloc)
cluster_method2.plot_cluster_weight(save_loc = saveloc)
saveloc = "./static/outputfiles/plots/" + resultId + "_11_choose_cluster_weight_GMM.png"
st.writePicture(saveloc)
outputText = cluster_method2.choose_clusters(thr_GMM) # thr_GMM = 0.2
st.write(outputText)

st.wf3step5()

saveloc = "./static/outputfiles/plots/" + resultId + "_12_conformation_KM.png"
st.writePicture(saveloc) # TO DO
ensemble_location1 = "./static/outputfiles/" + resultId + "_res_ensemble_kmeans"
outputText = cluster_method1.choose_conformations()
st.write(outputText)
cluster_method1.extract_conformations(ensemble_location1)

saveloc = "./static/outputfiles/plots/" + resultId + "_13_conformation_GMM.png"
st.writePicture(saveloc) # TO DO
ensemble_location2 = "./static/outputfiles/" + resultId + "_res_ensemble_gmm"
outputText = cluster_method2.choose_conformations()
st.write(outputText)
outputText = cluster_method2.extract_conformations(ensemble_location2)
st.write(outputText)

print("WF3 has been finished.")

cmd = "python sendEmail.py \"" + email + "\" \"" + resultId + "\" \"" + title + "\" \"" +  str(upload_time) + "\"" 
print(cmd)

if send_email_flag:
    Popen(cmd, shell=True)

file1 = open("./templates/output_results/" + resultId + "_outputHTML.html","w")
file1.write(st.HTML)
file1.close()


print(f"Analysis {resultId} has been finished.")
