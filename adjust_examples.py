GLPK_DIR="path/to/GLPK"
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
examples = ["combi_example", "combi_example_faults", "gene_distributed"]
for example in examples:
    pfilein = open(str(dir_path)+ "/distributedcombigrid/examples/" + example + "/Makefile" ,'r')
    temp = pfilein.read()
    pfilein.close()
    temp = temp.replace('$(SGPP)', str(dir_path))
    temp = temp.replace('$(GLPK)', GLPK_DIR)
    pfileout = open(str(dir_path)+ "/distributedcombigrid/examples/" + example + "/Makefile" ,'w')
    pfileout.write(temp)
    pfileout.close()

