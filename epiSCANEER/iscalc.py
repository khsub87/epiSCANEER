import sys, os, shutil
from Bio import AlignIO

from dhpylib import coe

# Global Function: get percentile rank
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
        
# Process of CN-Calculation        
class ProcCN:
    def __init__(self, proc_msa, coe_out_path, query_id, cn_cutoff=2.0, coe_algorithm="McBASC"):
        # Options (derivated from ProcMsa)
        self.proc_msa = proc_msa 
        self.aln_file_path = proc_msa.output_aln
        self.alignment = proc_msa.parse()
        self.coe_out_path = coe_out_path
        if not os.path.isabs(coe_out_path):
            self.coe_out_path = os.path.join(os.path.abspath('.'), coe_out_path)
        
        # Options (coe and cn options)
        self.query_id = query_id
        self.cn_cutoff = cn_cutoff
        self.coe_algorithm = coe_algorithm
        
        # Output
        self.CNScores = None
        self.CNPercentiles = None
        
    def calc(self):
        # Get Percentile Value of CNScores    
        def percentile_cn(self, raw_data, res_pos_dic, query_len):
            cn_dic = {}
            num_cutoff = int(self.cn_cutoff * query_len)
            
            raw_data.sort(key=lambda item:item[2], reverse=True)
            for (_ri, _rj, _) in raw_data[:num_cutoff]:
                ri, rj = (res_pos_dic[_ri], res_pos_dic[_rj])
                cn_dic[ri] = cn_dic.get(ri, 0) + 1
                cn_dic[rj] = cn_dic.get(rj, 0) + 1
                
            rList, rDic = get_percentile_rank(cn_dic.values())
            result = {}
            for a in cn_dic:
                result[a] = rDic[cn_dic[a]]
            return (cn_dic, result)                
        
        # calculation CN
        coe_in_path = self.aln_file_path + '_cn'
        error_path = self.aln_file_path + '_error'
        res_pos_dic, query_len = self.proc_msa.makeCNinput(coe_in_path, self.query_id, cutoff=0.2)
        Coe = coe.ProcCoe(coe_in_path, self.coe_out_path, error_path,self.coe_algorithm)
        Coe.run()
        cn_raw_data = Coe.parse()
        Coe.convertResPos(res_pos_dic)
        self.CNScores, self.CNPercentiles = percentile_cn(self, cn_raw_data, res_pos_dic, query_len)
        self.query_len = query_len
        return (self.CNScores, self.CNPercentiles)
        
# Process of conservation score calculation
class ProcCS:
    def __init__(self, proc_msa, output_res, query_id):
        # Options (derivated from ProcMsa)
        self.proc_msa = proc_msa
        self.aln_file_path = proc_msa.output_aln
        self.tree_file_path = proc_msa.output_tree
        self.alignment = proc_msa.parse()
        
        # Options (output)
        self.res_file_path = output_res
        self.query_id = query_id
        self.CSScores = None
        self.CSPercentiles = None
        
    def calc(self):
        def percentile_cs(raw_data):
            pr_con_data = {}
            aList = []
            for aData in raw_data:
                aList.append(aData[1])
            rList, rDic = get_percentile_rank(aList, bReverse=True)
            for aData in raw_data:
                pr_con_data[aData[0]] = rDic[aData[1]]
            return pr_con_data
        
        # Single-computer based Calculation Strategy           
        #else:
        options = {}
        options['fast_mode'] = True
        options['res_output'] = self.res_file_path
        options['tree_input'] = None
        options['query_id'] = self.query_id
        options['ori_output'] = self.res_file_path + '_ori'
        options['tree_output'] = None
        
        cs_in_path = self.aln_file_path + '_cs'
        self.proc_msa.makeCSinput(cs_in_path)
        res = rate4site.ProcRes(cs_in_path, options)
        res.run()
        cs_raw_data = res.parse()
        
        self.CSScores, self.CSPercentiles = (cs_raw_data, percentile_cs(cs_raw_data))
        return (self.CSScores, self.CSPercentiles)

# Process of Quantitive Integration Score calculation
class ProcIS:
    def __init__(self, proc_cn, proc_cs=None):
        self.proc_cn = proc_cn 
        self.proc_cs = proc_cs
        self.result = []
        
    def calc(self):
        CNScores, CNPercentiles, query_len = (self.proc_cn.CNScores, self.proc_cn.CNPercentiles, self.proc_cn.query_len)
        if self.proc_cs != None:
           CSScores, CSPercentiles = (self.proc_cs.CSScores, self.proc_cs.CSPercentiles)
        for res in range(0, query_len):
            if self.proc_cs != None:
               ISScore = CSPercentiles[res] * CNPercentiles.get(res, 0)
               CSScore = CSScores[res][1]
               CSPercentile = CSPercentiles[res]
            else:
               ISScore = 'N/A'
               CSScore = 'N/A'
               CSPercentile = 'N/A'

            self.result.append([str(res + 1),\
                                str(ISScore),\
                                str(CSScore),\
                                str(CSPercentile),\
                                str(CNScores.get(res, 0)),\
                                str(CNPercentiles.get(res, 0))])
        #self.result.sort(key=lambda item:item[0], reverse=False)
    
    def write(self, is_output):
        f = open(is_output, 'w')
        print >> f, '\t'.join(['res', 'IS', 'CS', 'CS_P', 'CN', 'CN_P'])
        for cell in self.result:
            print >> f, '\t'.join(cell)
        f.close()
