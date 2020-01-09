#!/usr/bin/env python

class TreeNode(object):
    def __init__(self):
        self.Value=None
        self.Type = None
        self.Filename = None
        self.isWritten = False
        self.Leaves=[]
        self.isTraversed=False

    def output(self):
        print "%s %s: %s"%(self.Type,self.Value,", ".join(self.Leaves))

    def appendLeaf(self,leaf):
        self.Leaves.append(leaf)

    def hasLeaves(self):
        if len(self.Leaves)>0:
            return True
        else:
            return False

    def hasSubtree(self):
        if len(self.Leaves)>0:
            if type(self.Leaves[0])==TreeNode:
                return True
        return False
    
    def setWritten(self):
        self.isWritten=True

    def isYetWritten(self):
        return self.isWritten

    def isModule(self):
        return self.Type=="M"

    def isProgram(self):
        return self.Type=="P"

def parseFortranSource(filename):
    result=[]

    srcfile=open(filename,"rt")
    for line in srcfile:
        #print line,
            
        mo = re.match(r"^\s*module\s+(\w+)(\s*$|\s*,\s*only)",line,re.I)
        if mo:
            # we found a module definition
            result.append(TreeNode())
            result[-1].Value=mo.group(1).lower()
            result[-1].Filename=filename
            result[-1].Type = "M"
            #print "Module definition found: ", result[-1].Value
            continue

        mo = re.match(r"^\s*end\s+module\s+",line,re.I)
        if mo:
            # module definition end
            continue

        mo =re.match(r"^\s*use\s+(\w+)",line,re.I)
        if mo:
            # we found a dependency
            #print "use statement found: ",mo.group(1),result[-1]
            result[-1].appendLeaf(mo.group(1).lower())
            continue

        mo = re.match(r"^\s*program\s+(\w+)",line,re.I)
        if mo:
            # main program found
            result.append(TreeNode())
            result[-1].Value=mo.group(1).lower()
            result[-1].Filename=filename
            result[-1].Type="P"
            continue

        mo = re.match(r"^\s*subroutine\s+(\w+)",line,re.I)
        if mo:
            # main program found
            result.append(TreeNode())
            result[-1].Value=mo.group(1).lower()
            result[-1].Filename=filename
            result[-1].Type="S"
            continue

        mo = re.match(r"^\s*\w*\s*function\s+(\w+)",line,re.I)
        if mo:
            # main program found
            result.append(TreeNode())
            result[-1].Value=mo.group(1).lower()
            result[-1].Filename=filename
            result[-1].Type="F"
            continue
            
    srcfile.close()
    return result

def findNode(name,nodelist):
    for node in nodelist:
        if (node.Value==name):
            return node
    return None

def buildSubtree(thisnode,nodes):
    # build a subtree 
    if not thisnode.hasSubtree():
        print "\nBEGIN subtree of node ",thisnode.Value
        leaflist = []
        for ileaf in range(len(thisnode.Leaves)):
            print "Getting %u. leaf from %u leaves."%(ileaf,len(thisnode.Leaves))
            #print thisnode.Leaves
            leaf = thisnode.Leaves.pop(0)
            #print "leaf = ",leaf,"\nnodes = ",[x.Value for x in nodes]
            node = findNode(leaf,nodes)
            if not node:
                # no module does exist ==> external
                #print "Adding an external node (Value=%s) to nodes"%leaf
                if with_external_nodes:
                    external_node=TreeNode()
                    external_node.Value=leaf
                    external_node.Filename="external_filename(%s)"%leaf
                    nodes.append(external_node)
                    leaflist.append(external_node)
            else:
                buildSubtree(node,nodes)
                leaflist.append(node)
        thisnode.Leaves=leaflist
        print "END subtree %s"%thisnode.Value
    else:
        print "Subtree of node %s already exists."%thisnode.Value

def writeNodesToFile(thisnode,dotfile):
    towrite=[]
    for subtree in thisnode.Leaves:
        if not(subtree.isYetWritten()):
            towrite.append(subtree.Value)
            subtree.setWritten()

    if len(towrite)>1:
        dotfile.write('{ rank=same; ')
        for string in towrite:
            dotfile.write('%s; '%string)
        dotfile.write(' }\n')

    for subtree in thisnode.Leaves:
        if not(subtree.isTraversed):
            #print "traversing subtree = ",subtree.Value
            writeNodesToFile(subtree,dotfile)
            subtree.isTraversed=True

def writeToFile(thisnode,dotfile):
    print "Writing subtrees for node ",thisnode.Value
    for subtree in thisnode.Leaves:
        dotfile.write('\t%s -> %s;\n'%(thisnode.Value,subtree.Value))

    for subtree in thisnode.Leaves:
        if not(subtree.isTraversed):
            writeToFile(subtree,dotfile)
            subtree.isTraversed=True

def resetAllNodes(nodes):
    for node in nodes:
        node.isWritten=False
        node.isTraversed=False

def writeMakeDependencies(thisnode,makedepfile):
    print "Writing make dependencies for node ",thisnode.Value
    for subtree in thisnode.Leaves:
        # construct the makefile line by line
        # target object file should reside in $(OBJDIR)
        # src is in $(SRCDIR)
        if not (subtree.Filename == thisnode.Filename):
            target = "$(OBJDIR)/%s.o"%os.path.splitext(os.path.basename(thisnode.Filename))[0]
            dependency = "$(OBJDIR)/%s.o"%os.path.splitext(os.path.basename(subtree.Filename))[0]
            makedepfile.write('%s: %s\n'%(target,dependency))

    for subtree in thisnode.Leaves:
        if not(subtree.isTraversed):
            writeMakeDependencies(subtree,makedepfile)
            subtree.isTraversed=True

def getLinkLine(thisnode,makedepfile):
    if not thisnode.isYetWritten():
        towrite=["$(OBJDIR)/%s.o"%os.path.splitext(os.path.basename(thisnode.Filename))[0]]
        thisnode.setWritten()
    else:
        towrite=[]

    for subtree in thisnode.Leaves:
        if not(subtree.isYetWritten()):
            towrite.append("$(OBJDIR)/%s.o"%os.path.splitext(os.path.basename(subtree.Filename))[0])
            subtree.setWritten()

    #if len(towrite)>1:
    #    for string in towrite:
    #        makedepfile.write('%s '%string)

    for subtree in thisnode.Leaves:
        if not(subtree.isTraversed):
            #print "traversing subtree = ",subtree.Value
            towrite.extend(getLinkLine(subtree,makedepfile))
            subtree.isTraversed=True

    return towrite

import sys
import re
import os.path
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-r","--root",dest="rootnode",default="gene")
parser.add_option("-e","--with-external-nodes",dest="with_external_nodes",
                  action="store_true",default=False)

(options,args) = parser.parse_args()

filelist = args
rootname=options.rootnode
with_external_nodes = options.with_external_nodes

nodes = []
for thisfile in filelist:
    # is the suffix a Fortran90 suffix?
    if re.search(r"\.(f90|F90|F03|f03|F95|f95)\s*$",os.path.splitext(thisfile)[1]):
        print "parsing ",thisfile
        resultlist = parseFortranSource(thisfile)
        #print "resultlist = ",[x.Value for x in resultlist]
        for res in resultlist:
            if res.isModule() or res.isProgram():
                nodes.append(res)

dotfile=open('use_deps.dot',"w")
dotfile.write("digraph G {\n")
for node in nodes:
    for leaf in node.Leaves:
        dotfile.write("\t%s -> %s;\n"%(node.Value,leaf))

dotfile.write("}\n")
dotfile.close()

# we now have all modules, programs, subroutines, function which
# contains use statements in a list 
#rootname = 'df_nonlinear_term_mod'
#print "Looking for rootname %s\n"%rootname
#print nodes
for node in nodes:
    if node.Value==rootname:
        rootnode = node
        break
        #tree.setRoot(node)
        
# now we have the root of the tree fixed
#print rootnode.Leaves
print [x.Value for x in nodes]

buildSubtree(rootnode,nodes)

# starting with rootnode, we have the tree structure
dotfile=open('use_tree.dot',"w")
dotfile.write("digraph G {\n")
dotfile.write('orientation=landscape;\n')
dotfile.write('size="15.5,11";\n')
dotfile.write('ratio=fill;\n')
dotfile.write('margin="0.5,0.5";\n')
writeNodesToFile(rootnode,dotfile)
resetAllNodes(nodes)
writeToFile(rootnode,dotfile)
dotfile.write("}\n")
dotfile.close()

resetAllNodes(nodes)
# now we also write a dependency file in makefile style
#print rootnode
#print [x.Filename for x in rootnode.Leaves]
makedepfile = open('dependencies.mk',"w")
linkline = "%s.o"%os.path.splitext(os.path.basename(rootnode.Filename))[0]
linkline = ' \\\n'.join(["\t%s"%x for x in getLinkLine(rootnode,makedepfile)])
makedepfile.write("LINK_OBJ = %s\n\n"%linkline)
#for linkobj in linkline:
#    makedepfile.write("\t%s \\\n"%linkobj)

resetAllNodes(nodes)
writeMakeDependencies(rootnode,makedepfile)
makedepfile.close()
print linkline

