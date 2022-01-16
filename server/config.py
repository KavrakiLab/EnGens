#uploader 
# ------------------------------------------------------------------------------------------ #


from pickle import FALSE, TRUE

from numpy import true_divide


domain_url = "http://127.0.0.1:5000" #domain for the URL
upload_folder = 'static/uploads/'   #UPLOAD FOLDER PATH

email_regex_string= r'([A-Za-z0-9]+[.-_])*[A-Za-z0-9]+@[A-Za-z0-9-]+(\.[A-Z|a-z]{2,})+'   
ALLOWED_EXTENSIONS = set(['pdb', 'xtc', 'dcd', 'xtx'])
secret_key = "engens-123456789"

# user = 'root'
# password = 'admin12345678'
# host = '127.0.0.1'
# database = 'kavraki'

#workflow1
# ------------------------------------------------------------------------------------------ #

# 1 - load full trajectory from the files and visualize it; 
# 2 - load a subset of trajectory using a list of atom indices;
# 3 - load a subset of trajectory using the atom selection string
engen_object_config = 3

# load a subset of trajectory using a list of atom indices;
first_N_atoms = 200

# known residues of the binding site to start from
binding_site_residues = [6, 56, 60, 63, 64, 78, 81, 82, 85, 86, 89, 94, 150, 164, 165, 166, 171, 174, 178, 179, 182, 183, 243, 246, 247, 249, 250, 253, 261, 262, 263, 264, 267, 268, 271, 274, 275]


#workflow2
# ------------------------------------------------------------------------------------------ #

reduction_type = "TICA" #can be "TICA"|"HDE"|"PCA"
chosen_lag = 10

var_thr = 70
comp_num = 16


#workflow3
# ------------------------------------------------------------------------------------------ #
clustering = "KM" # options: KM | GMM

clusterRange = range(2, 10)

# KM
param_index_KM = 2
param_index_GMM = 0

thr_KM = 0.2
thr_GMM = 0.2


# emails
# ------------------------------------------------------------------------------------------ #

send_email_flag = TRUE
#send_email_flag = FALSE

# password = "KavrakiLab123!"
# port = 465
# smtp_server = "smtp.gmail.com"
# sender_email = "automated.message.from.the.lab@gmail.com"  # Enter your address