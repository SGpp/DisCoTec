#!/usr/bin/env python

import AutoparFile
import sys
import numpy as np
import matplotlib.colors
import matplotlib.cm as cm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.axisartist as axisartist

show_bar_graph=False
# first read in the autopar.dat given as argument

af=AutoparFile.AutoparFile()
af.read(sys.argv[1])

if (len(sys.argv)>2):
    af2=AutoparFile.AutoparFile()
    af2.read(sys.argv[2])
    #print af2.data_array

#print af.data_array

(best_par,best_pv)=af.getBestParallelization()
print "Best parallelization: ",best_par,", with perf_vec = ",best_pv

# now we have to convert the data to a numpy array
if (show_bar_graph):
    normalize = matplotlib.colors.Normalize()
    normalize.autoscale_None(af.mem_array.flatten())
    scmap = cm.ScalarMappable(norm=normalize, cmap=cm.jet)
    color_array =  scmap.to_rgba(af.mem_array.flatten())

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.bar3d(af.parall_array.flatten(),af.perfvec_array.flatten(),np.zeros_like(data_array).flatten(),0.5,0.5,af.data_array.flatten(),color=color_array)

    plt.show()


######################## plot the data in an array plot #########
if (len(sys.argv)==2):
    # print data_array
    plt.matshow(np.transpose(af.data_array))
    # Labeling of the axes
    plt.xticks(range(af.paral_index),af.paral_label,rotation='vertical')

    # initialize the list with None values, because we need the whole
    # list for nonsequential access with the next loop
    pv_label=[]
    for x in range(af.perf_vec_index):
        pv_label.append(None)
        
    for key in af.perf_vec_label:
        pv_label[af.perf_vec_label[key]]=key

    plt.yticks(range(af.perf_vec_index),pv_label)

    # add a colorbar
    plt.colorbar()

    fig2=plt.gcf()
    fig2.savefig('autopar_plot.pdf',bbox_inches='tight')
    # plt.show()

########### plot two plots aligned ##########
if (len(sys.argv)>2):
    # we have the files af and af2
    # first get crop the files to the same minimal length
    # assuming that they only differ at the end
    len1=len(af.data_array)
    len2=len(af2.data_array)
    #print len1,len2
    minlen=min(len1,len2)
    cutarray1=np.transpose(af.data_array[0:minlen,:])
    cutarray2=np.transpose(af2.data_array[0:minlen,:])
    cut_paral_label=af.paral_label[0:minlen]

    # normalize both to the minimum runtime
    min1 = np.nanmin(cutarray1)
    min2 = np.nanmin(cutarray2)
    print "minvals = ", min1,min2
    normarray1 = cutarray1/min1
    normarray2 = cutarray2/min2

    max1 = np.nanmax(normarray1)
    max2 = np.nanmax(normarray2)
    print "maxvals after normalization = ", max1, max2
    #plt.subplots(2,1,sharex=True)
    # Calculate aspect ratio of the data

    # initialize the pv_label list with None values, because we need the whole
    # list for nonsequential access with the next loop
    up_pv_label=[]
    for x in range(af.perf_vec_index):
        up_pv_label.append(None)
        
    for key in af.perf_vec_label:
        up_pv_label[af.perf_vec_label[key]]=key

    lo_pv_label=[]
    for x in range(af2.perf_vec_index):
        lo_pv_label.append(None)
        
    for key in af2.perf_vec_label:
        lo_pv_label[af2.perf_vec_label[key]]=key

    
    data_aspect_ratio1=float(cutarray1.shape[1])/cutarray1.shape[0]
    data_aspect_ratio2=float(cutarray2.shape[1])/cutarray2.shape[0]
    #print data_aspect_ratio1,data_aspect_ratio2

    width = 0.8
    up_height = width/data_aspect_ratio1
    lo_height = width/data_aspect_ratio2
    between = 0.02
    v_center = 0.25
    left = 0.1
    up_bottom = v_center+0.5*between
    lo_bottom = v_center-0.5*between-lo_height
    co_bottom = 0.05
    co_height = 0.05

    fig=plt.figure(1)
    upax=fig.add_axes([left,up_bottom,width,up_height])
    upimg = upax.imshow(normarray1,aspect='auto',interpolation='nearest')
    upax.xaxis.set_ticks(range(minlen))
    upax.xaxis.set_ticklabels(cut_paral_label,size='small',rotation='vertical')
    for tick in upax.xaxis.get_major_ticks():
        tick.tick1On=False
        tick.tick2On=True
        tick.label1On=False
        tick.label2On=True
        tick.label2.set(rotation='vertical',size='x-small')
    #upax.xaxis.set_label('Different parallelizations')
    #upax.xaxis.set_label_coords(0.5,0.5,transform=upax.transAxes)
    #upax.xaxis.label.set(position=(0.5,0.5),visible=True)

    for tick in upax.yaxis.get_major_ticks():
        tick.label1On=False
        tick.label2On=False
    upax.set_ylabel('Linux cluster',rotation='horizontal')

    #upax.yaxis.set_ticks(range(af.perf_vec_index))
    #upax.yaxis.set_ticklabels(up_pv_label,size='x-small')
    #print "upax = ", upax.get_position()
    #print upax.get_aspect()

    loax=fig.add_axes([left,lo_bottom,width,lo_height])
    loimg = loax.imshow(normarray2,aspect='auto',interpolation='nearest')
    for tick in loax.xaxis.get_major_ticks():
        tick.tick1On=False
        tick.tick2On=False
        tick.label1On=False
        tick.label2On=False
     
    for tick in loax.yaxis.get_major_ticks():
        tick.label1On=False
        tick.label2On=False

    loax.set_ylabel('Power 6',rotation='horizontal')
    print loax.yaxis.get_label()

    #loax.yaxis.set_ticks(range(af2.perf_vec_index))
    #loax.yaxis.set_ticklabels(lo_pv_label,size='x-small')
    #print "lowax = ", lowax.get_position()
    #print lowax.get_aspect()
    coax=fig.add_axes([left,co_bottom,width,co_height])
    fig.colorbar(upimg,cax=coax,orientation='horizontal')

    fig.savefig('double_autopar.pdf', bbox_inches='tight')

    print fig.axes
    #plt.show()

    
