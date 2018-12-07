#!/usr/bin/env python
import os 
import sys
import os.path
CC = sys.argv[1]
FC = sys.argv[2]
cc = sys.argv[3]


dir_path = os.path.dirname(os.path.realpath(__file__))
GLPK_DIR= str(dir_path) + "/glpk"
examples = ["combi_example", "combi_example_faults", "gene_distributed", "gene_distributed_linear"]
for example in examples:
    pfilein = open(str(dir_path)+ "/distributedcombigrid/examples/" + example + "/Makefile.template" ,'r')
    temp = pfilein.read()
    pfilein.close()
    temp = temp.replace('$(SGPP)', str(dir_path))
    temp = temp.replace('$(GLPK)', GLPK_DIR)
    temp = temp.replace('$(CC)', CC)
    temp = temp.replace('$(FC)', FC)
    temp = temp.replace('$(cc)', cc)
    pfileout = open(str(dir_path)+ "/distributedcombigrid/examples/" + example + "/Makefile" ,'w')
    pfileout.write(temp)
    pfileout.close()
gene_versions = ["gene_mgr", "gene-non-linear", "gene-dev-exahd"]
for geneVersion in gene_versions:
    fnameList = [str(dir_path) + "/"+ geneVersion + "/makefiles/new_machine/new_machine", str(dir_path) + "/"+ geneVersion + "/makefiles/hazel_hen/hazel_hen"]
    for fname in fnameList:
        if(os.path.isfile(fname + ".template")):
            with open(fname + ".template",'r') as pfilein:
                temp = pfilein.read()
            temp = temp.replace('$(SGPP)', str(dir_path))
            temp = temp.replace('$(GLPK)', GLPK_DIR)
            temp = temp.replace('$(CC)', CC)
            temp = temp.replace('$(FC)', FC)
            temp = temp.replace('$(cc)', cc)
            with open(fname + ".mk" ,'w') as pfileout:
                pfileout.write(temp)
pfilein = open(str(dir_path)+ "/distributedcombigrid/examples/gene_distributed/preproc.py" ,'r')
temp = pfilein.read()
pfilein.close()
temp = temp.replace('$(SGPP)', str(dir_path))
pfileout = open(str(dir_path)+ "/distributedcombigrid/examples/gene_distributed/preproc.py" ,'w')
pfileout.write(temp)
pfileout.close()

pfilein = open(str(dir_path)+ "/distributedcombigrid/examples/gene_distributed_linear/preproc.py" ,'r')
temp = pfilein.read()
pfilein.close()
temp = temp.replace('$(SGPP)', str(dir_path))
pfileout = open(str(dir_path)+ "/distributedcombigrid/examples/gene_distributed_linear/preproc.py" ,'w')
pfileout.write(temp)
pfileout.close()

