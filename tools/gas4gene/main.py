"""This module provide a main function that reads data
from files specified in the argument line."""

# external packages and modules
import sys
import argparse
import numpy as np


#--------------------------------------------------------------------
try:
    import matplotlib.pyplot as plt
# ploting of the surface of the distribution function
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
# set flag
    plotting_loaded = True
except ImportError:
    plotting_loaded = False
    print "no matplotlib packages! no plotting!"
#--------------------------------------------------------------------


#local packages and modules
import readers
from gridvx import Gridvx
from gridwx import Gridwx
from profiles import ProfilesData

def main():
    """generate adaptive grid for GENE."""

    blk_adp_grid = True

    # parsing command line arguments
    parser = argparse.ArgumentParser(description = "generate adaptive static grid for GENE")
    parser.add_argument("parameters_file_name", help = "name of parameters file")
    parser.add_argument("--nvxblks", type = int, \
                        help = "number of blocks for v||, x block-structured subgrid")
    parser.add_argument("--nwxblks", type = int, \
                        help = "number of blocks for mu, x block-structured subgrid")
    erf_group = parser.add_mutually_exclusive_group()
    erf_group.add_argument("--erft", nargs = 2, help = \
                           "read erf table file name and choose precision level" \
                           + "(x, v|| subgrid)")
    erf_group.add_argument("--probvx", type = float, help = \
                           "compute confidence interval for v||, x subgrid constraction")
    erf_group.add_argument("--probwx", type = float, help = \
                           "compute confidence interval for mu, x subgrid constraction")
    args = parser.parse_args()

    # get information about profiles
    parameters_str = readers.get_config_str(args.parameters_file_name)
    profile_dict = readers.get_profile_info(parameters_str)

    if (profile_dict['prof_type'] == -1) :
        # read profile with temperatures and densities
        profile_array = readers.read_profile(profile_dict['name'])
    else:
        profile_array = None

    # get precision and confidence interval saved in prec_intv_pair_vx
    # for x, v|| subgrid
    # first element [0] is tolerance or precision for truncated distribution
    # second element [1] is n from confidence itnerval expression (-n*sigma, n*sigma)
    if args.erft:
        prec_intv_pair_vx = readers.get_essential_interval_normal_dist(
            erf_table_file_name = args.erft[0], plevel = int(args.erft[1]))
    elif args.probvx:
        prec_intv_pair_vx = readers.get_essential_interval_normal_dist(
            tolerance = args.probvx)
    elif ('prob_vx' in profile_dict) and (profile_dict['prob_vx'] != None):
        prec_intv_pair_vx = \
            readers.get_essential_interval_normal_dist(
            tolerance = profile_dict['prob_vx'])
    else:
        print "--erft or --prob_vx flag should be specified or information", \
              "about probability should be written in the parameters file"
        return 1

    # get precision and confidence interval saved in prec_intv_pair_wx
    # for x, w subgrid
    # first element [0] is tolerance or precision for truncated distribution
    # second element [1] is scaled by Tprf/Tx0 interval for subgrid constracution
    if args.probwx:
        prec_intv_pair_wx = readers.get_essential_interval_exp_dist(
            tolerance = args.probwx)
    elif ('prob_wx' in profile_dict) and (profile_dict['prob_wx'] != None):
        prec_intv_pair_wx = readers.get_essential_interval_exp_dist(
            tolerance = profile_dict['prob_wx'])
    else :
        print "--prob_wx flag should be specified or information", \
              "about probability should be written in the parameters file"
        return 1

    #number of blocks for v||, x block-structured subgrid generation
    if args.nvxblks:
        num_of_blks_vx = args.nvxblks
    elif ('n_vx_blks' in profile_dict) and (profile_dict['n_vx_blks'] != None):
        num_of_blks_vx = profile_dict['n_vx_blks']
    else:
        print "--nvxblks should be spceciefed or n_vx_blks should", \
              "be written in the parameters file"
        return 1

    #number of blocks for mu, x block-structured subgrid generation
    if args.nwxblks:
        num_of_blks_wx = args.nwxblks
    elif ('n_wx_blks' in profile_dict) and (profile_dict['n_wx_blks'] != None):
        num_of_blks_wx = profile_dict['n_wx_blks']
    else:
        print "--nwxblks should be spceciefed or n_wx_blks should", \
              "be written in the parameters file"
        return 1

    # read the method name for finding an optimal x, v|| block-grid structure
    if ('opt_mthd_vx' in profile_dict) and (profile_dict['opt_mthd_vx'] != None):
        opt_mthd_vx = profile_dict['opt_mthd_vx']
    else:
        print "default optimization method is chosen for the x, v|| subgrid,", \
              "to change the method specify its name in the parameters file"
        opt_mthd_vx = 'sum'

    # read the method name for finding an optimal x, mu block-grid structure
    if ('opt_mthd_wx' in profile_dict) and (profile_dict['opt_mthd_wx'] != None):
        opt_mthd_wx = profile_dict['opt_mthd_wx']
    else:
        print "default optimization method is chosen for the x, v|| subgrid,", \
              "to change the method specify its name in the parameters file"
        opt_mthd_wx = 'sum'
    
    # initializing the profiles from the data available
    profilesDT = ProfilesData(profile_dict, profile_array)
        
    # initialize v||, x subgrid object
    blk_grid_vx = Gridvx(profilesDT, prec_intv_pair_vx)

    # initialize mu, x subgrid object
    blk_grid_wx = Gridwx(profilesDT, prec_intv_pair_wx)

    # init optimization methods dictionary (might be extended)
    methods_dict_vx = { 'sum' : blk_grid_vx.getOptBGridDiscreteStr, \
                        'integral' : blk_grid_vx.getOptBGridSmoothStr }

    methods_dict_wx = { 'sum' : blk_grid_wx.getOptBGridDiscreteStr, \
                        'integral' : blk_grid_wx.getOptBGridSmoothStr }

    blk_mks_r_v = methods_dict_vx[opt_mthd_vx](num_of_blks_vx)
    blk_mks_v = blk_grid_vx.getVelocityBlkMks(blk_mks_r_v)

    #blk_mks_r_w = methods_dict_wx[opt_mthd_wx](num_of_blks_wx)
    blk_mks_r_w = blk_mks_r_v
    blk_mks_w = blk_grid_wx.getMagneticMomentBlkMks(blk_mks_r_w)
    
    # optimization of the grid structure (merging blocks)
    blk_grid_wx.optimizeBlkStr(blk_mks_r_w, blk_mks_w)
    blk_mks_r_w    = blk_grid_wx.blk_mks_r_w
    blk_mks_w      = blk_grid_wx.blk_mks_w
    num_of_blks_wx = blk_grid_wx.n_wx_blks

    #save vx block-structured grid results into the parameter file
    parameters_str = \
        readers.add_vx_blk_structure_info(parameters_str, num_of_blks_vx, \
                                             prec_intv_pair_vx[0],
                                             opt_mthd_vx,
                                             blk_mks_r_v, blk_mks_v)

    # readers.save_config_str(args.parameters_file_name, \
    #                             readers.get_output_config_str(parameters_str))
    
    #save wx block-structured grid results into the parameter file
    parameters_str = \
        readers.add_wx_blk_structure_info(parameters_str, num_of_blks_wx, \
                                             prec_intv_pair_wx[0],
                                             opt_mthd_wx,
                                             blk_mks_r_w, blk_mks_w)

    readers.save_config_str(args.parameters_file_name, \
                                readers.get_output_config_str(parameters_str))
    
    #------------------------------ visualization part ----------------------------------

    if(plotting_loaded):

        # visualization of temperature and density profiles
        # plot_temp_dens_profiles(profilesDT)

        #------------------------ visualization of grids --------------------------------

        # get coordinate lines ticks for all blocks
        if(blk_adp_grid == True):
            (r_v_blks, v_blks) = blk_grid_vx.getBlkTicks_badptv(blk_mks_r_v, blk_mks_v)
        else:
            (r_v_blks, v_blks) = blk_grid_vx.getBlkTicks(blk_mks_r_v, blk_mks_v)
        (r_w_blks, w_blks) = blk_grid_wx.getBlkTicks()
        
        plot_vx_wx_subgrids(blk_grid_vx, blk_grid_wx, prec_intv_pair_vx, prec_intv_pair_wx,
                            blk_mks_r_v, blk_mks_v, blk_mks_r_w, blk_mks_w,
                            r_v_blks, v_blks, r_w_blks, w_blks, profilesDT)

        # unify v-x and w-x block-structured subgrids
        (r_vwx_blks, v_vwx_blks, w_vwx_blks) = \
                     constructVWXGrid(r_v_blks, v_blks, r_w_blks, w_blks)

        plot_vwx_subgrid(r_vwx_blks, v_vwx_blks, w_vwx_blks)
        
        # find number of points in the constant resolution grid
        if(blk_adp_grid == True):
            (np_rg, np_vx, np_wx, np_vwx) = estimate_number_of_points_badptv( \
                profilesDT, r_v_blks, v_blks, r_w_blks, w_blks,\
                r_vwx_blks, v_vwx_blks, w_vwx_blks)
        else:
            (np_rg, np_vx, np_wx, np_vwx) = estimate_number_of_points( \
                profilesDT, r_v_blks, v_blks, r_w_blks, w_blks,\
                r_vwx_blks, v_vwx_blks, w_vwx_blks)
            
        print "number of blocks ->"
        print "in v|| direction    : ", len(r_v_blks)
        print "in mu direction     : ", len(r_w_blks)
        print "in combinned version: ", len(r_vwx_blks)
        print "number of points in the regular grid: ", int(np_rg)
        print "number of points in the vx b-s grid : ", int(np_vx),  " -> ", np_vx/np_rg
        print "number of points in the wx b-s grid : ", int(np_wx),  " -> ", np_wx/np_rg
        print "number of points in the vwx b-s grid: ", int(np_vwx), " -> ", np_vwx/np_rg
                
        # plot w-v-x subgrid region boundary
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        w = np.arange(0, blk_mks_w[0], blk_mks_w[0]/20.0)
        w = np.hstack((w, blk_mks_w[0]))
        v = np.arange(-blk_mks_v[0], blk_mks_v[0], blk_mks_v[0]/20.0)
        v = np.hstack((v, blk_mks_v[0]))
        w, v = np.meshgrid(w, v)

        # find temperature for a given v range
        sigma = np.divide(v,prec_intv_pair_vx[1])
        tv  = np.divide(np.multiply(sigma,sigma),2)
        # find temperature for a given w range
        tw  = np.divide(w,prec_intv_pair_wx[1])

        # general value
        tvw = np.maximum(tv,tw)
        # scale temperature
        tvw = np.multiply(tvw, profilesDT.ref_temperature)
        # find flux surface label (radial distance)
        x = profilesDT.getInvTemperature(tvw)

        surf = ax.plot_surface(w, v, x, rstride=1, cstride=1, cmap=cm.jet, \
                                   linewidth=0, antialiased=False)
        
        ax.set_zlim(profilesDT.xval_a_min, profilesDT.xval_a_max)
    
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
        #   fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set_xlabel('mu')
        ax.set_ylabel('v||')
        ax.set_zlabel('x')
        
        plt.show()
             
        # ploting surface of distribution function

        # fig = plt.figure()
        # ax = fig.gca(projection='3d')
        # X = np.arange(velocity_range[0]/2, velocity_range[-1]/2, 0.05)
        # Y = np.arange(profilesDT.xval_a_min, profilesDT.xval_a_max, 0.05)
        # X, Y = np.meshgrid(X, Y)
        # Z = blk_grid_vx.getDistrSurf(X, Y)
        # surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet, \
        #                            linewidth=0, antialiased=False)
        # ax.set_zlim(0, 1.5)
    
        # ax.zaxis.set_major_locator(LinearLocator(10))
        # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
        # #    fig.colorbar(surf, shrink=0.5, aspect=5)
        # ax.set_xlabel('parallel velocity')
        # ax.set_ylabel('flux surface label')
        # ax.set_zlabel('pdf')
                    
        # fig_exp_contr = plt.figure()
        # ax = fig_exp_contr.gca(projection='3d')
        # X = np.arange(0, blk_mks_w[0]/2, 0.05)
        # Y = np.arange(profilesDT.xval_a_min, profilesDT.xval_a_max, 0.05)
        # X, Y = np.meshgrid(X, Y)
        # Z = blk_grid_wx.getDistrSurf(X, Y)
        # surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet, \
        #                            linewidth=0, antialiased=False)
        # ax.set_zlim(0, 2.5)
    
        # ax.zaxis.set_major_locator(LinearLocator(10))
        # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
        # #    fig.colorbar(surf, shrink=0.5, aspect=5)
        # ax.set_xlabel('parallel velocity')
        # ax.set_ylabel('flux surface label')
        # ax.set_zlabel('pdf')

        # plt.figure(4)
        # plt.hold(True)

        # m = 10
        # n = 20

        # (x, y) = blk_grid_vx.getGridTfI(m, n)

        # plt.plot(velocity_range,radius_points,'go')

        # plt.title('parallel velocity component - radius phase subspace area \n' \
        #               + 'with ' + str(prec_intv_pair_vx[0]*100) + "% of all particles")

        # for i in range(0, m):
        #     plt.plot(x[i,:], y[i,:], 'b-')

        # for j in range(0, n):
        #     plt.plot(x[:,j], y[:,j], 'b-')
                    
        # plt.xlabel('parallel velocity component')
        # plt.ylabel('radial distance')

        # plt.hold(False)

        # plt.show()
        
    return 0

def constructVWXGrid(r_v_blks, v_blks, r_w_blks, w_blks):
    
    """combine separate structures for vx and wx subgrids to obtain vwx subgrid."""

    r_vwr_blks = []
    v_vwr_blks = []
    w_vwr_blks = []

    vb = True
    wb = True

    iv  = 0
    iw  = 0
    nw  = 0
    nv  = 0
    nvw = 0
    
    while ( vb or wb):
        if( nv+len(r_v_blks[iv]) < nw+len(r_w_blks[iw]) ):

            # print "iv: ", iv
            # print "iw: ", iw

            r_vwr_blk = r_v_blks[iv][nvw-nv:]
            r_vwr_blks.append(r_vwr_blk)
            v_vwr_blks.append(v_blks[iv])
            w_vwr_blks.append(w_blks[iw])
            
            nv += len(r_v_blks[iv])-1
            nvw = nv
            if(iv < len(r_v_blks)-1):
                iv += 1
            else:
                vb = False
        elif( nv+len(r_v_blks[iv]) > nw+len(r_w_blks[iw]) ):

            # print "iv: ", iv
            # print "iw: ", iw

            r_vw_blk = r_w_blks[iw][nvw-nw:]
            r_vwr_blks.append(r_vw_blk)
            v_vwr_blks.append(v_blks[iv])
            w_vwr_blks.append(w_blks[iw])
            
            nw += len(r_w_blks[iw])-1
            nvw = nw
            if(iw < len(r_w_blks)-1):
                iw += 1
            else:
                wb = False
        else:

            # print "iv: ", iv
            # print "iw: ", iw

            r_vw_blk = r_w_blks[iw][nvw-nw:]
            r_vwr_blks.append(r_vw_blk)
            v_vwr_blks.append(v_blks[iv])
            w_vwr_blks.append(w_blks[iw])
                        
            nw += len(r_w_blks[iw])-1
            nv += len(r_v_blks[iv])-1
            nvw = nw
            if(iw < len(r_w_blks)-1):
                iw += 1
            else:
                wb = False
            if(iv < len(r_v_blks)-1):
                iv += 1
            else:
                vb = False
                
    return (r_vwr_blks, v_vwr_blks, w_vwr_blks)

def plot_temp_dens_profiles(profilesDT):
    
    profile_array = profilesDT.profile_array
        
    plt.figure(figsize=(8, 10))

    # temperature profile

    plt.subplot(211)
    
    x = np.arange(profilesDT.xval_a_min, profilesDT.xval_a_max, .01)
    p = profilesDT.getTemperature(x)
    
    plt.plot(profile_array[:,0], profile_array[:,2], 'ro', x, p, 'k-')
    plt.title('temperature profile')
    plt.xlabel('flux surface label')
    plt.ylabel('temperature')
    
    # electrons concentration profile
    plt.subplot(212)
    
    x = np.arange(profilesDT.xval_a_min, profilesDT.xval_a_max, .01)
    d = profilesDT.getDensity(x)
    
    plt.plot(profile_array[:,1], profile_array[:,3], 'bs', x, d, 'b-')
    plt.title('electrons concentration profile')
    plt.xlabel('flux surface label')
    plt.ylabel('electrons concentration')

def plot_vx_wx_subgrids(blk_grid_vx, blk_grid_wx,
                        prec_intv_pair_vx, prec_intv_pair_wx,
                        blk_mks_r_v, blk_mks_v,
                        blk_mks_r_w, blk_mks_w,
                        r_v_blks, v_blks,
                        r_w_blks, w_blks,
                        profilesDT):

    profile_array = profilesDT.profile_array
    
    # determine parallel velocity range
    velocity_range = np.multiply(-prec_intv_pair_vx[1], blk_grid_vx.sigma)
    velocity_range = \
                   np.hstack((velocity_range, np.multiply(prec_intv_pair_vx[1], \
                                                          blk_grid_vx.sigma)[::-1]))
    # determine flux surface label (radius)
    radius_points = np.hstack((profile_array[:,0],profile_array[::-1,0]))
    
    plt.figure()
        
    plt.hold(True)
        
    bgrid_struct = blk_grid_vx.getBGridVizStructure(blk_mks_r_v)

    plt.plot(velocity_range,radius_points,'go',\
                 bgrid_struct[1,:], bgrid_struct[0,:], 'b-')
    # plot coordinate lines
    for ib in range(0, len(r_v_blks)):
        r_v_blk = r_v_blks[ib]
        v_blk = v_blks[ib]

        # horizontal lines
        n = r_v_blk.size
        for ih in range(0, n):
            plt.plot([v_blk[0],v_blk[-1]], [r_v_blk[ih],r_v_blk[ih]], 'y')

        # vertical lines
        n = v_blk.size
        for iv in range(0, n):
            plt.plot([v_blk[iv],v_blk[iv]], [r_v_blk[0],r_v_blk[-1]], 'y')

    plt.title('parallel velocity component - radius phase subspace area \n' \
              + 'with ' + str(prec_intv_pair_vx[0]*100) + "% of all particles")
    plt.xlabel('parallel velocity component')
    plt.ylabel('radial distance')

    plt.hold(False)

    plt.figure()
    plt.hold(True)

    bgrid_struct = blk_grid_wx.getBGridVizStructure(
        blk_mks_r_w, blk_mks_w)
        
    plt.plot(blk_grid_wx.mu, profile_array[:,0],'go',
             bgrid_struct[1,:], bgrid_struct[0,:], 'b-')

    # plot coordinate lines
    for ib in range(0, len(r_w_blks)):
        r_w_blk = r_w_blks[ib]
        w_blk = w_blks[ib]

        # horizontal lines
        n = r_w_blk.size
        for ih in range(0, n):
            plt.plot([w_blk[0],w_blk[-1]], [r_w_blk[ih],r_w_blk[ih]], 'y')
            
        # vertical lines
        n = w_blk.size
        for iv in range(0, n):
            plt.plot([w_blk[iv],w_blk[iv]], [r_w_blk[0],r_w_blk[-1]], 'y')

    plt.title('magnetic moment component - radius phase subspace area \n' \
              + 'with ' + str(prec_intv_pair_wx[0]*100) + "% of all particles")
    plt.xlabel('magnetic moment component')
    plt.ylabel('radial distance')

    plt.hold(False)

def plot_vwx_subgrid(r_vwx_blks, v_vwx_blks, w_vwx_blks):
            
     # ploting block-structured grid in x, v||, mu subspace
        
     fig = plt.figure()
     plt.hold(True)
        
     ax = fig.gca(projection='3d')
     #ax.w_xaxis.set_scale('log')
        
     # plot coordinate lines
     for ib in range(0, len(r_vwx_blks)):
         r_vwx_blk = r_vwx_blks[ib]
         v_vwx_blk = v_vwx_blks[ib]
         w_vwx_blk = w_vwx_blks[ib]
         n_r = r_vwx_blk.size
         n_v = v_vwx_blk.size
         n_w = w_vwx_blk.size

         # horizontal lines
         # just surface of the block!
            
         for iv in range(0, n_v):
             line = ax.plot(
                 [w_vwx_blk[0] ,w_vwx_blk[-1]],
                 [v_vwx_blk[iv],v_vwx_blk[iv]],
                 [r_vwx_blk[0],r_vwx_blk[0]],
                 'y')
             line = ax.plot(
                 [w_vwx_blk[0] ,w_vwx_blk[-1]],
                 [v_vwx_blk[iv],v_vwx_blk[iv]],
                 [r_vwx_blk[-1],r_vwx_blk[-1]],
                 'y')

         for iw in range(0, n_w):
             line = ax.plot(
                 [w_vwx_blk[iw],w_vwx_blk[iw]],
                 [v_vwx_blk[0] ,v_vwx_blk[-1]],
                 [r_vwx_blk[0],r_vwx_blk[0]],
                 'y')
             line = ax.plot(
                 [w_vwx_blk[iw],w_vwx_blk[iw]],
                 [v_vwx_blk[0] ,v_vwx_blk[-1]],
                 [r_vwx_blk[-1],r_vwx_blk[-1]],
                 'y')
            
         for ir in range(0, n_r):
             line = ax.plot(
                 [w_vwx_blk[0] ,w_vwx_blk[0]],
                 [v_vwx_blk[0] ,v_vwx_blk[-1]],
                 [r_vwx_blk[ir],r_vwx_blk[ir]],
                 'y')
             line = ax.plot(
                 [w_vwx_blk[-1],w_vwx_blk[-1]],
                 [v_vwx_blk[0] ,v_vwx_blk[-1]],
                 [r_vwx_blk[ir],r_vwx_blk[ir]],
                 'y')

         for ir in range(0, n_r):
             line = ax.plot(
                 [w_vwx_blk[0] ,w_vwx_blk[-1]],
                 [v_vwx_blk[0] ,v_vwx_blk[0]],
                 [r_vwx_blk[ir],r_vwx_blk[ir]],
                 'y')
             line = ax.plot(
                 [w_vwx_blk[0] ,w_vwx_blk[-1]],
                 [v_vwx_blk[-1],v_vwx_blk[-1]],
                 [r_vwx_blk[ir],r_vwx_blk[ir]],
                 'y')
              
         # vertical lines

         for iv in range(0, n_v):
             line = ax.plot(
                 [w_vwx_blk[0] ,w_vwx_blk[0]],
                 [v_vwx_blk[iv],v_vwx_blk[iv]],
                 [r_vwx_blk[0] ,r_vwx_blk[-1]],
                 'y')
             line = ax.plot(
                 [w_vwx_blk[-1],w_vwx_blk[-1]],
                 [v_vwx_blk[iv],v_vwx_blk[iv]],
                 [r_vwx_blk[0] ,r_vwx_blk[-1]],
                 'y')

         for iw in range(0, n_w):
             line = ax.plot(
                 [w_vwx_blk[iw],w_vwx_blk[iw]],
                 [v_vwx_blk[0] ,v_vwx_blk[0]],
                 [r_vwx_blk[0] ,r_vwx_blk[-1]],
                 'y')
             line = ax.plot(
                 [w_vwx_blk[iw],w_vwx_blk[iw]],
                 [v_vwx_blk[-1],v_vwx_blk[-1]],
                 [r_vwx_blk[0] ,r_vwx_blk[-1]],
                 'y')
                
     ax.set_xlabel('mu')
     ax.set_ylabel('v||')
     ax.set_zlabel('x')


def estimate_number_of_points(\
    profilesDT, r_v_blks, v_blks, r_w_blks, w_blks, \
    r_vwx_blks, v_vwx_blks, w_vwx_blks):
    # full regular grid
    np_rg = profilesDT.nx0*profilesDT.nky0*profilesDT.nz0\
            *profilesDT.nv0*profilesDT.nw0*profilesDT.n_spec
    
    # vx block-structured grid
    n_vx_blks = len(r_v_blks)
    np_vx = r_v_blks[0].size*v_blks[0].size
    for i in range(1,n_vx_blks):
        np_vx += (r_v_blks[i].size-1)*v_blks[i].size
    np_vx *= profilesDT.nky0*profilesDT.nz0*profilesDT.nw0*profilesDT.n_spec

    # wx block-structured grid
    n_wx_blks = len(r_w_blks)
    np_wx = r_w_blks[0].size*w_blks[0].size
    for i in range(1,n_wx_blks):
        np_wx += (r_w_blks[i].size-1)*w_blks[i].size
    np_wx *= profilesDT.nky0*profilesDT.nz0*profilesDT.nv0*profilesDT.n_spec

    # vwx block-structured grid
    n_vwx_blks = len(r_vwx_blks)
    np_vwx = r_w_blks[0].size*v_blks[0].size*w_blks[0].size
    for i in range(1,n_vwx_blks):
        np_vwx += (r_vwx_blks[i].size-1)*v_vwx_blks[i].size*w_vwx_blks[i].size
    np_vwx *= profilesDT.nz0*profilesDT.n_spec

    return (np_rg, np_vx, np_wx, np_vwx)

def estimate_number_of_points_badptv(\
    profilesDT, r_v_blks, v_blks, r_w_blks, w_blks,\
    r_vwx_blks, v_vwx_blks, w_vwx_blks):
    # full regular grid
    np_rg = profilesDT.nx0*profilesDT.nky0*profilesDT.nz0\
            *profilesDT.nv0*profilesDT.nw0*profilesDT.n_spec
    
    # vx block-structured grid
    n_vx_blks = len(r_v_blks)
    np_vx = r_v_blks[0].size*v_blks[0].size
    for i in range(1,n_vx_blks):
        np_vx += r_v_blks[i].size*v_blks[i].size
    np_vx *= profilesDT.nky0*profilesDT.nz0*profilesDT.nw0*profilesDT.n_spec

    # wx block-structured grid
    n_wx_blks = len(r_w_blks)
    np_wx = r_w_blks[0].size*w_blks[0].size
    for i in range(1,n_wx_blks):
        np_wx += r_w_blks[i].size*w_blks[i].size
    np_wx *= profilesDT.nky0*profilesDT.nz0*profilesDT.nv0*profilesDT.n_spec

    # vwx block-structured grid
    n_vwx_blks = len(r_vwx_blks)
    np_vwx = r_w_blks[0].size*v_blks[0].size*w_blks[0].size
    for i in range(1,n_vwx_blks):
        np_vwx += r_vwx_blks[i].size*v_vwx_blks[i].size*w_vwx_blks[i].size
    np_vwx *= profilesDT.nz0*profilesDT.n_spec

    return (np_rg, np_vx, np_wx, np_vwx)
       
if __name__ == "__main__":
    sys.exit(main())
