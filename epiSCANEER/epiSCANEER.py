# IS Calculation Script
import os, glob, sys, math, shutil
from dhpylib import msa, coe, blast, iscalc
import numpy as np
import networkx as nx
import csv

def build_coevolution_score_dic(base_path, prefix, cn_kind):
    output_dic = {}
    f = open(os.path.join(base_path, "%s.coe_out_%s" %(prefix, cn_kind)))
    f.readline()
    for line in f.readlines():
        line = line.strip().split("\t")
        output_dic[(line[0], line[1])] = float(line[2])
    f.close()
    return output_dic

def get_percentile_rank(aList, bReverse=False):
    step = 1.0 / (len(aList) - 1)
    bList = []
    for v in aList:
        bList.append(v)
    bList.sort(reverse=bReverse)
    rDic, cur = ({}, 0)
    for v in bList:
        if v not in rDic:
            rDic[v] = cur 
        cur += step
    rList = []
    for v in aList:
        rList.append(rDic[v])
    return (rList, rDic)

def calc_CN(sorted_key, len_seq, len_threshold):
    CN_dic = {}
    for res in range(1, len_seq+1):
        CN_dic[res] = 0 
    for i in range(int(len_seq*len_threshold)):
        res1, res2 = sorted_key[i]
        CN_dic[int(res1)] += 1
        CN_dic[int(res2)] += 1
    values = CN_dic.values()
    tmp_list, per_dic = get_percentile_rank(values)
    per_CN_dic = {}
    for res in CN_dic.keys():
        per_CN_dic[res] = per_dic[CN_dic[res]]
    return CN_dic, per_CN_dic

def build_coupling_dic(sorted_key, len_seq, len_threshold):
    coupling_dic = {}
    for res in range(1, len_seq+1):
        coupling_dic[res] = []
    for i in range(int(len_seq*len_threshold)):
        res1, res2 = sorted_key[i]
        coupling_dic[int(res1)].append(int(res2))
        coupling_dic[int(res2)].append(int(res1))
    return coupling_dic

def build_msa(msa_path, prefix):
	tmp_dic = {}
	f = open(msa_path)
	for line in f.readlines():
		line = line.strip().split()
		if len(line) == 2:
			if prefix[:30] in line[0]: line[0] = prefix
			if not line[0] in tmp_dic:
				tmp_dic[line[0]] = ""
			tmp_dic[line[0]] += line[1]
	f.close()
	
	msa_dic = {}
	gene = prefix
	for i in range(len(tmp_dic[gene])):
		if tmp_dic[gene][i] != "-":
			for key in tmp_dic.keys():
				if not key in msa_dic:
					msa_dic[key] = ""
				try: msa_dic[key] += tmp_dic[key][i]
				except IndexError: msa_dic[key] += "-"
	return msa_dic

def get_co_evolving_graph(coenet_path):
    with open(coenet_path) as f:
        f.readline()
        G = nx.Graph()
        for line in f:
            line = line.strip().split("\t")
            G.add_edge(int(line[0]), int(line[1]))
    return G

def get_frequency(msa_dic, pos1, pos2, AA1, AA2):
    residue_list = [
        seq[pos1-1] + seq[pos2-1]
        for seq in msa_dic.values()
        if "-" not in (seq[pos1-1] + seq[pos2-1])
    ]
    return residue_list.count(AA1 + AA2) / float(len(residue_list))


def get_num(msa_dic, pos1, pos2, AA1, AA2):
    residue_list = [
        seq[pos1-1] + seq[pos2-1]
        for seq in msa_dic.values()
        if "-" not in (seq[pos1-1] + seq[pos2-1])
    ]
    return residue_list.count(AA1 + AA2)


def get_multiSCI(msa_dic, seq, G, pos1, pos2, AA1, AA2):
    edges = set()
    for pos in [pos1, pos2]:
        for neighbor in G.neighbors(pos):
            edges.add(tuple(sorted([pos, neighbor])))

    SCI_list = []
    for node1, node2 in edges:
        if node1 in (pos1, pos2):
            AA_node1 = AA1 if node1 == pos1 else AA2
        else:
            AA_node1 = seq[node1-1]

        if node2 in (pos1, pos2):
            AA_node2 = AA1 if node2 == pos1 else AA2
        else:
            AA_node2 = seq[node2-1]

        WT_num = get_num(msa_dic, node1, node2, seq[node1-1], seq[node2-1])
        MUT_num = get_num(msa_dic, node1, node2, AA_node1, AA_node2)

        if seq[node1-1] == AA_node1 and seq[node2-1] == AA_node2:
            SCI_list.append(0)
        else:
            SCI_list.append(math.log((MUT_num + 1) / float(WT_num)))

    SCI_list.extend(np.zeros(G.number_of_edges() - len(SCI_list)))
    return np.mean(SCI_list)


def get_AA_variance(msa_dic, pos):
    return {
        seq[pos-1]
        for seq in msa_dic.values()
        if seq[pos-1] != "-"
    }


def final_SCI(prefix):
    infile = os.path.join(base_pth, f"{prefix}.doubleSCI_mp")
    outfile = os.path.join(base_pth, f"{prefix}.epiSCI_final")

    rows = []

    with open(infile, newline="") as f:
        reader = csv.DictReader(f, delimiter='\t')
        reader.fieldnames = [h.strip() for h in reader.fieldnames]
        for row in reader:
            row["MUT_SCI"] = float(row["MUT_SCI"])
            row["Double"] = row["Double"].lower() == "true"
            row["pos1"] = int(row["pos1"])
            row["pos2"] = int(row["pos2"])
            rows.append(row)

    sci_vals = [r["MUT_SCI"] for r in rows]
    sci_min = min(sci_vals)
    sci_max = max(sci_vals)

    for r in rows:
        r["SCI_norm"] = (r["MUT_SCI"] - sci_min) / (sci_max - sci_min)

    best_per_pair = {}

    for r in rows:
        if not r["Double"]:
            continue

        key = (r["pos1"], r["pos2"])
        if key not in best_per_pair or r["SCI_norm"] > best_per_pair[key]["SCI_norm"]:
            best_per_pair[key] = r

    final_rows = list(best_per_pair.values())

    fieldnames = final_rows[0].keys()
    with open(outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames,delimiter="\t")
        writer.writeheader()
        writer.writerows(final_rows)

def calc_epiSCI(base_pth, msa_dic, pm, prefix):
	# Get CN (default: McBASC)
	print ("CS/CN calculation of %s" % prefix)
	coe_output_mcbasc = os.path.join(base_pth, "%s.coe_out_mcbasc" % prefix)
	pcn = iscalc.ProcCN(pm, coe_output_mcbasc, pm.result[0].id, cn_cutoff=2.0, coe_algorithm="McBASCCovariance")
	cn_result = pcn.calc()

	tmp_coe = build_coevolution_score_dic(base_pth, prefix, 'mcbasc')
	sorted_coe = sorted(tmp_coe, key=lambda x: tmp_coe[x], reverse=True)
	CN_dic, CN_per_dic = calc_CN(sorted_coe, pcn.query_len, 2.0)
	
	# Print CN results
	CNScores, query_len = (pcn.CNScores, pcn.query_len)
	f = open(os.path.join(base_pth, "%s.cn") %(prefix), 'w')
	print( "\t".join(['res', 'CN_McBasc', 'CN_per']),file=f)
	for res in range(1, query_len+1):
		output = [res, CN_dic.get(res,0), CN_per_dic.get(res, 0)]
		print("\t".join(map(str,output)),file=f)
	f.close()

	# Build and extract coevolutionary network
	coenet_dic = build_coupling_dic(sorted_coe, pcn.query_len, 2.0)
	f = open(os.path.join(base_pth, "%s.coenet") %(prefix), 'w')
	print("res1\tres2",file=f)
	for res1 in coenet_dic.keys():
		for res2 in coenet_dic[res1]:
			if res1 > res2:
				print("\t".join(map(str,[res1,res2])),file=f)
	f.close()

	#calc_Double_mutation
	out_path = os.path.join(base_pth, "%s.doubleSCI_mp") %(prefix)

	msa_dic = build_msa(os.path.join(base_pth, "%s.aln" %prefix), prefix)
	seq=msa_dic[prefix]

	G = get_co_evolving_graph(os.path.join(base_pth, "%s.coenet") %(prefix))
	edges = list(G.edges())

	print ("gene epiSCI value of %s..." % prefix)
	with open(out_path, "w") as fo:
		header = [
			"pos1", "pos2", "WT_AA1", "WT_AA2",
			"WT_freq", "MUT_AA1", "MUT_AA2",
			"MUT_SCI", "Double"
		]
		print("\t".join(header), file=fo)

		for num, (res1, res2) in enumerate(edges):
			AA_set1 = get_AA_variance(msa_dic, res1)
			AA_set2 = get_AA_variance(msa_dic, res2)

			AA_A, AA_B = seq[res1 - 1], seq[res2 - 1]
			freq = get_frequency(msa_dic, res1, res2, AA_A, AA_B)

			for AA1 in AA_set1:
				for AA2 in AA_set2:
					multiSCI = get_multiSCI(
						msa_dic, seq, G, res1, res2, AA1, AA2
					)
					row = [
						res1, res2, AA_A, AA_B,
						freq, AA1, AA2,
						multiSCI,
						AA_A != AA1 and AA_B != AA2
					]
					print("\t".join(map(str, row)), file=fo)
	
	final_SCI(prefix)

