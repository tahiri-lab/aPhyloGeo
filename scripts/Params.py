import yaml
from yaml.loader import SafeLoader


# We open the params.yaml file and put it in the params variable
with open('./scripts/params.yaml') as f:
    params = yaml.load(f, Loader=SafeLoader)

bootstrap_threshold =  params["bootstrap_threshold"]
ls_threshold =         params["ls_threshold"]
window_size =          params["window_size"]
step_size =            params["step_size"]
data_names =           ['ALLSKY_SFC_SW_DWN_newick', 
                                'T2M_newick', 
                                'QV2M_newick', 
                                'PRECTOTCORR_newick', 
                                'WS10M_newick']
reference_gene_file =  params["reference_gene_dir"]+'/'+params["reference_gene_file"]
reference_gene_filename = params["reference_gene_file"]

file_name =            params["file_name"]
specimen =             params["specimen"]
names      =           params["names"]
makeDebugFiles=        params["makeDebugFiles"]
bootstrapAmount=        params["bootstrap_amount"]