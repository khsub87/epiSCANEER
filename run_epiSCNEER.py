from epiSCANEER.epiSCANEER import *
from epiSCANEER import msa, coe, iscalc

input_path = "./input/"
output_path = "./output/"
if not os.path.isdir(output_path):
	os.makedirs(output_path)

file_list = glob.glob(os.path.join(input_path, "*.aln"))
for msa_path in file_list:
	msa_file = os.path.basename(msa_path)
	prefix = msa_file.split('.')[0]

	base_pth = os.path.join(output_path, prefix)
	if not os.path.isdir(base_pth):
		os.makedirs(base_pth)
	elif len(os.listdir(base_pth)) == 5:
		print "already calculated! skip!"
		continue

	# Loading MSA
	pm = msa.ProcMsa("tmp", msa_path, "tmp", "tmp")
	pm.parse()
	msa_dic = build_msa(msa_path)

	# Construct co-evolutionary network and Calculating SCI
    calc_epiSCI(base_pth, msa_dic, pm, prefix)