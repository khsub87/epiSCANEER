import os, sys
from dhpylib import lib_path

class ProcCoe:
    # input aln, out coe file path are required
    def __init__(self, aln_file, coe_file, error_file,algorithm="McBASC"):
        self.algorithm = algorithm
        self.input_aln_file = os.path.abspath(aln_file)
        self.output_coe_file = os.path.abspath(coe_file)
        self.result = []
    
    # Available Algorithm: McBASC, ELSC, SCA, MI, OMES
    def run(self):
        if sys.platform == 'win32': # Windows Environment
            jexePath = 'java.exe'
            cmd_apd = ' & '
        else: # Linux Environment
            jexePath = os.path.join(lib_path, "jre/linux/bin/java")
            jexePath = "/usr/bin/java"
            cmd_apd = ';'
        jClassPath = os.path.join(lib_path, "coe")
        jClassAlgorithm = "covariance.algorithms.%s" % self.algorithm
        exejClass = "%s %s %s %s" % (jexePath, jClassAlgorithm, self.input_aln_file, self.output_coe_file) 
        os.system(cmd_apd.join(['cd %s' % jClassPath, exejClass]))
        
    # Parse coe-calculation result
    def parse(self):
        f = open(self.output_coe_file, 'r')
        #f.next()
        f.readline()
        for line in f.readlines():
            fields = line.split()
            self.result.append((int(fields[0]), int(fields[1]), float(fields[2])))
        f.close()
        return self.result
    
    def convertResPos(self, res_pos_dic):
        f = open(self.output_coe_file, 'w')
        print ("residue_1\tresidue_2\tscore",file=f)
        for cell in self.result:
            print ('%d\t%d\t%s' %\
             (res_pos_dic[cell[0]]+1, res_pos_dic[cell[1]]+1, str(cell[2])), file=f)
        f.close()
            
