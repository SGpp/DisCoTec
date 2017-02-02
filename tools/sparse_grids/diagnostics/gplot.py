import os
from math import *
from numpy import *

class graph:
    def __init__(self,lineformat,xval,yval,zval=[]):
        self.format=lineformat
        self.x=xval
        self.y=yval
        self.z=zval    

    def write(self,pfile):

        if ((type(self.x)!=list) and (type(self.x)!=ndarray)):
            #1D
            line='%15g   %15g\n' % (self.x,self.y)
            pfile.write(line)
            pfile.write('end\n')
        elif ((type(self.x[0])!=list) and (type(self.x[0])!=ndarray)):
            #2D
            if self.z==[]:
                #no color coding
                for ind in range(len(self.x)):
                    line='%15g   %15g\n' % (self.x[ind],self.y[ind])
                    pfile.write(line)
            else:
                #with color coding
                for ind in range(len(self.x)):
                    line='%15g   %15g  %15g\n' % (self.x[ind],self.y[ind],self.z[ind])
                    pfile.write(line)
            pfile.write('end\n')
        else:
            #3D
            for ind1 in range(len(self.x)):
                for ind2 in range(len(self.x[ind1])):
                    line='%15g   %15g  %15g\n' % (self.x[ind1][ind2],self.y[ind1][ind2],self.z[ind1][ind2])
                    pfile.write(line)
                pfile.write('\n')
            pfile.write('end\n')

class plot:
    def __init__(self,obj_list,xlabel,ylabel,zlabel='',title='',loglog=False,xrange=[],yrange=[],contonly=False):
        self.obj_list=obj_list
        self.xlabel=xlabel
        self.ylabel=ylabel
        self.zlabel=zlabel
        self.title=title
        self.loglog=loglog
        self.xrange=xrange
        self.yrange=yrange
        self.contonly=contonly
        
    def write(self,pfile):
        line="set title '%s'\n" % self.title
        pfile.write(line)
        line="set xlabel '%s'\n" % self.xlabel
        pfile.write(line)
        line="set ylabel '%s'\n" % self.ylabel
        pfile.write(line)
        if self.xrange==[]:
            line='set autoscale x\n'
            pfile.write(line)
        else:
            line='set xrange [%f:%f]\n' % (self.xrange[0],self.xrange[1])
            pfile.write(line)
        if self.yrange==[]:
            line='set autoscale y\n'
            pfile.write(line)
        else:
            line='set yrange [%f:%f]\n' % (self.yrange[0],self.yrange[1])
            pfile.write(line)

        if self.loglog:
            pfile.write('set log xy\n') 
        if self.zlabel =='':
            #2D plot
            line='plot '
            for i in range(len(self.obj_list)):
                if i>0:
                    line+=', '
                if self.obj_list[i].z==[]:
                    line +="'-' using 1:2 %s" % self.obj_list[i].format
                else:
                    line +="'-' using 1:2:3 %s" % self.obj_list[i].format
        else:
            #3D plot
            if self.contonly:
                line+='set pm3d map\n'
            else:
                line="set zlabel '%s'\n" % self.zlabel
                line+='set pm3d\n'
            pfile.write(line)
            line='splot '
            for i in range(len(self.obj_list)):
                if i>0:
                    line+=', '
                line +="'-' using 1:2:3 %s" % self.obj_list[i].format

        line +='\n'
        pfile.write(line)
        for obj in self.obj_list:
            obj.write(pfile)

class multiplot:
    def __init__(self,plotlist,title='',pause=-1):
        self.plotlist=plotlist
        self.title=title
        self.pause=pause
        self.nplots=len(self.plotlist)
        self.nrows=int(ceil(sqrt(self.nplots)))
        self.nlines=int(ceil(float(self.nplots)/self.nrows))
        self.scro=1.0/self.nrows
        self.scli=1.0/self.nlines

    def write(self,pfile):
        pfile.write('set multiplot\n')
        pfile.write('set size %f,%f\n' % (self.scro,self.scli))
        for i in range(self.nplots):
            pfile.write('set origin %f,%f\n' % (i%self.nrows*self.scro,(self.nlines-1-i/self.nrows)*self.scli))
            self.plotlist[i].write(pfile)
                    
        pfile.write('set nomultiplot\n')
        pfile.write('pause %d\n' % self.pause)
            
def write_plotfile(mplotlist,title):
    pfile=open('%s.plt' % title,'w')
    pfile.write('# Gnuplot script file for plotting gene/sg++ results\n')
    for mplot in mplotlist:
        mplot.write(pfile)
    
    pfile.close

def show(mplotlist,title):
    write_plotfile(mplotlist,title)
    if mplotlist[0].nplots==1:
         os.system('gnuplot %s.plt' % title)
    else:
        os.system('gnuplot -geometry 1200x900 %s.plt' % title)

