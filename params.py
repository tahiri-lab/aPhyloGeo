#file_name = 'The_37_climate.csv'
file_name = '5seq/geo.csv'

specimen = 'id'   #"Please enter the name of the colum containing the specimens names: "

names = ['id', 'ALLSKY_SFC_SW_DWN', 'T2M', 'PRECTOTCORR', 'QV2M', 'WS10M']

bootstrap_threshold = 0
rf_threshold = 100
window_size = 200
step_size = 100
data_names = ['ALLSKY_SFC_SW_DWN_newick', 'T2M_newick', 'QV2M_newick', 'PRECTOTCORR_newick', 'WS10M_newick']

#reference_gene_file = 'datasets/The_37seq.fasta'

reference_gene_file = 'datasets/5seq/seq.fasta'
