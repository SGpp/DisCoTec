#!/usr/bin/env python

import sys
import re
import os.path
from optparse import OptionParser


def search_for_output(startdir,jobid):
    filename_to_find = "GENE.%u.out"%jobid
    for dir in os.walk(startdir):
        if re.search(r"prob\d+",dir[0]):
            for file in dir[2]:
                if file==filename_to_find:
                    return os.path.join(dir[0],filename_to_find)
                #endif
            #endfor
        #endif
    #endfor
#end def

parser = OptionParser("extract_perf.py [Optionen] 'GENE Output Dateien mit perf Teil'")
parser.add_option("-s","--sort",dest="SORTKEY",
                  help="Sortkey for sort operation, use 'calls', 'inc/exc_time/percent/MFlops'",
                  default="exc_time")
parser.add_option("-c","--compare",dest="COMPLABELS",
                  help="colon separated list of perf labels or the special label 'datasize', which will be compared in a multiple run comparison",
                  default="t_loop")
parser.add_option("-m","--metric",dest="METRIC",
                  help="metric for the comparison of multiple runs",
                  default="exc_time")
parser.add_option("-t","--context", dest="CONTEXT",
                  help="name of the interesting context, if nothing given, the global context is used",
                  default="global context")

(options,args) = parser.parse_args()

if len(args) < 1:
    parser.error("At least one gene output filename has to be given as argument.")
#end if

# we have to interpret the options.COMPLABELS as this is a colon separated list
complabel_list = options.COMPLABELS.split(':')
input_context_name = options.CONTEXT
input_context_name.strip()

perflist = []
for gene_out_file in args:
    mo = re.search(r"^(\d+)$",gene_out_file)
    if mo:
        # the argument was just the jobid
        jobid = int(mo.group(1))
        # we have to go through all ../prob* directories to look for
        # GENE.jobid.out files
        gene_out_file = search_for_output("..",jobid)
        print "Working on ",gene_out_file
    #end if

    # Get the jobid and the path out of the passed filename
    #input_directory = os.path.dirname(gene_out_file)
    mo = re.search(r"GENE.(\d+).out",os.path.basename(gene_out_file))
    if mo:
        jobid = int(mo.group(1))
    #end if
    mo = re.search(r"geneerr.log",os.path.basename(gene_out_file))
    if mo:
        jobid = 0
    #end if


    infile = open(gene_out_file)
    extract = False
    for line in infile:
        mo = re.search(r"Choice for parallelization:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)",line)
        if mo:
            #outfile.write(line)
            perflist.append({"jobid":jobid,"n_procs_s":int(mo.group(1)),
                             "n_procs_v":int(mo.group(2)),
                             "n_procs_w":int(mo.group(3)),
                             "n_procs_x":int(mo.group(4)),
                             "n_procs_y":int(mo.group(5)),
                             "n_procs_z":int(mo.group(6)),
                             "nblocks":1})
            blocks = []
        #end if

        mo = re.search(r"^\s*nblocks\s*=\s*(\d+)",line)
        if mo:
            perflist[-1]['nblocks'] =int(mo.group(1))

        mo = re.search(r"^=== Context:\s+(.+)\s+===",line)
        if mo:
            context_name = mo.group(1)
            context_name.strip()
            if context_name==input_context_name:
                in_right_context = True
            else:
                in_right_context = False

        if re.search(r"\s*Inclusive\s*Exclusive",line) and in_right_context:
            extract = True
        #end if

        if extract:
            #outfile.write(line)
            #print line
            #mo = re.search(r"(\w{,8})\s*(\d+)\s+([0-9.]+)\s+([0-9.]+)",line)
            mo = re.search(r"(\w{,8})\s*(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)",line)
            if mo:
                #print mo.group(1),mo.group(2),mo.group(3)
                blocks.append({'label':mo.group(1),
                               'calls':int(mo.group(2)),
                               'inc_time':float(mo.group(3)),
                               'inc_percent':float(mo.group(4)),
                               'inc_MFlops':float(mo.group(5)),
                               'exc_time':float(mo.group(6)),
                               'exc_percent':float(mo.group(7)),
                               'exc_MFlops':float(mo.group(8))})
            #end if
        #end if

        mo = re.search(r"\s*Size of data segment used by the program:\s*([0-9.]+)\s*MB",line)
        if mo:
            extract = False
            perflist[-1]["blocks"]=blocks
            perflist[-1]["datasize"]=float(mo.group(1))
        #end if
    #end for

    #outfile.close()
    infile.close()
#end for

print "Number of analyzed runs is ",len(perflist)
    
# perflist is a list of all files given as input. Each file is a dictionary with the elements 
# jobid, all processes of the parallelization  and another dictionary blocks, which contains
# the different instrumented perf blocks
    
# now we implement the actions for the command line options
first_run = True
for this_run in perflist:
    # sorting
    #print "Sorting %u with sortkey %s"%(this_run["jobid"],options.SORTKEY)
    sorted_file = open("sorted.%u.perf"%this_run["jobid"],"w")
    blocks = this_run["blocks"]
    blocks.sort(key=(lambda a:a[options.SORTKEY]),reverse=True )
    sorted_file.write("sorted by %s\n\n"%options.SORTKEY)
    
    for block in blocks:
        sorted_file.write("%8s %10u %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n"
                          %(block["label"],block["calls"],
                            block["inc_time"],block["inc_percent"],block["inc_MFlops"],
                            block["exc_time"],block["exc_percent"],block["exc_MFlops"]))
    #end for
    sorted_file.close()
    
    # comparison for multiple runs
    #print "comparison of the blocks",complabel_list," for the metric",options.METRIC
    blocklabel_list = [x["label"] for x in blocks]
    if first_run:
        print "Metric printed is ",options.METRIC
        print "jobid   s  v  w  x  y  z nblocks",
        for label in complabel_list:
            print "%10s"%label,
        print
        first_run = False
    #end if
    print "%6s %2u %2u %2u %2u %2u %2u %5u "%(this_run["jobid"],this_run["n_procs_s"],this_run["n_procs_v"],
                                      this_run["n_procs_w"],this_run["n_procs_x"],
                                      this_run["n_procs_y"],this_run["n_procs_z"],this_run["nblocks"]),
    
    for complabel in complabel_list:
        if complabel=="datasize":
            print "%10g"%this_run["datasize"],
        else:
            try:
                comp_index = blocklabel_list.index(complabel)
            except ValueError:
                print "          ",
            else:
                print "%10g"%blocks[comp_index][options.METRIC],
            #end try
        #end if
    #end for
    print
#end for
