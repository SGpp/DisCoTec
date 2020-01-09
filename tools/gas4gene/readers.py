"""This modules provides functions to read and write profiles of temperatures and densities of TOKAMAK or STELLARATOR cross-section"""

import ConfigParser
import numpy as np
from scipy.special import erfinv
from math import sqrt, log
from string import replace
from io import BytesIO
from os import linesep

def read_erf_table(erf_table_file_name):

    """read error function table from a file"""
    
    erf_table_file = open(erf_table_file_name,'r')
    print erf_table_file

    # skipping comment line
    erf_table_file.readline()

    # loop over meaningfull data
    erf_table = np.array([float(val) for val in erf_table_file.readline().split()])
    for line in erf_table_file:
        erf_table = \
            np.vstack((erf_table, np.array([float(val) for val in line.split()])))

    erf_table_file.close()

    return erf_table

def get_essential_interval_normal_dist(
    erf_table_file_name = None, plevel = 0, tolerance = 0.98):

    """returns a tuple with precision and scaled by standard deviation confidence interval \
    for normal distribution"""

    if erf_table_file_name is None:
        return (tolerance, sqrt(2)*erfinv(tolerance))
    else:
        erf_table = read_erf_table(erf_table_file_name)
        return (erf_table[plevel, 0], erf_table[plevel, 1])

def get_essential_interval_exp_dist(tolerance = 0.98):

    """ returns a tuple with precision and interval scaled by 1/(profile temperature) \
    for exponential distribution"""

    return (tolerance, -log(1-tolerance))

def read_profile(profile_file_name):

    """ read temperature and dencities profile from a file"""
    
    profile_file = open(profile_file_name, 'r')
    print profile_file

    # skipping comment lines
    profile_file.readline()
    profile_file.readline()

    # loop over meaningful data
    profile_array = np.array([float(val) for val in profile_file.readline().split()])
    for line in profile_file:
        profile_array = \
            np.vstack((profile_array, np.array([float(val) for val in line.split()])))

    profile_file.close()
        
    return profile_array


def get_config_str(parameters_file_name):

    """read parameters file and convert the content into string which can be parsed with ConfigParser"""

    # read parameters file and save the content in the string
    parameters_file = open(parameters_file_name)
    parameters_str = parameters_file.read()
    parameters_file.close()

    # get sections from the parameters record
    sections = [word for word in parameters_str.split() if word.startswith('&')]

    # replace &section by [section], so it can be parsed with ConfigParser
    for section in sections:
        parameters_str = replace(parameters_str,section,'['+section[1:]+']')

    # remove all empty lines
    parameters_str = linesep.join([s for s in parameters_str.splitlines() if s])

    return parameters_str


def save_config_str(parameters_file_name, parameters_str):

   """save configuration string into a file """

   # open a parameters file and save the content in it
   parameters_file = open(parameters_file_name, 'w')
   parameters_file.write(parameters_str)
   parameters_file.close()
    

def get_output_config_str(parameters_str):

    """convert the parameters string modified by ConfigParser into a standard GENE parameters file"""

    # get sections from the parameters string
    parameters_out_str = parameters_str

    sections = [word for word in parameters_out_str.split() if word.startswith('[')]

    # replace [section] by &section, so it can be parsed with ConfigParser
    for section in sections:
        parameters_out_str = replace(parameters_out_str, section, '&'+section[1:-1])

    # add empty line between sections
    parameters_out_str = replace(parameters_out_str, '/\n', '/\n\n')

    return parameters_out_str    

def get_profile_type_pair(parameters_str):
    
    """parse parameters string and return the value of the profile type and experimental profile file name (meaningless for analytical profiles) in a tuple"""

    config = ConfigParser.RawConfigParser(allow_no_value = True)
    config.readfp(BytesIO(parameters_str))

    prof_type = config.getint('species', 'prof_type')
    name = config.get('species', 'name')
    name = name[1:-1]

    return (prof_type,'profile_' + name)


def get_profile_info(parameters_str):

    """parse parameters string and retrieve information necessary for grid generation"""

    config = ConfigParser.RawConfigParser(allow_no_value = True)
    config.readfp(BytesIO(parameters_str))

    rhostar = config.getfloat('geometry', 'rhostar')
    x0 = config.getfloat('box', 'x0')
    lx = config.getfloat('box', 'lx')
    nx0 = config.getfloat('box', 'nx0')
    
    main_dict = dict( \
        prof_type = config.getint('species', 'prof_type'), \
        name = 'profile_' + config.get('species', 'name')[1:-1], \
        kappa_T = config.getfloat('species', 'kappa_T'), \
        LT_center = config.getfloat('species', 'LT_center'), \
        LT_width = config.getfloat('species', 'LT_width'), \
        delta_x_T = config.getfloat('species', 'delta_x_T') if \
            config.has_option('species', 'delta_x_T') else None, \
        kappa_n = config.getfloat('species', 'kappa_n'), \
        Ln_center = config.getfloat('species', 'Ln_center'), \
        Ln_width = config.getfloat('species', 'Ln_width'), \
        delta_x_n = config.getfloat('species', 'delta_x_n') if \
            config.has_option('species', 'delta_x_n') else None, \
        minor_r = config.getfloat('geometry', 'minor_r'), \
        major_R = config.getfloat('geometry', 'major_R'), \
        xval_a_min = x0 - 0.5*rhostar*lx, \
        xval_a_max = x0 + 0.5*rhostar*lx, \
        dx = lx*rhostar/(nx0-1),\
        x0 = x0, \
        nx0 = nx0, \
        nky0 = config.getfloat('box', 'nky0'), \
        nz0 = config.getfloat('box', 'nz0'), \
        nv0 = config.getfloat('box', 'nv0'), \
        nw0 = config.getfloat('box', 'nw0'), \
        n_spec = config.getfloat('box', 'n_spec')
        )

    aux_dict = dict()
    if config.has_section('adaptivity'):
        aux_dict = dict( \
            n_vx_blks = config.getint('adaptivity', 'n_vx_blks') if \
                config.has_option('adaptivity', 'n_vx_blks') else None, \
            prob_vx = config.getfloat('adaptivity', 'prob_vx') if \
                config.has_option('adaptivity', 'prob_vx') else None, \
            n_wx_blks = config.getint('adaptivity', 'n_wx_blks') if \
                config.has_option('adaptivity', 'n_wx_blks') else None, \
            prob_wx = config.getfloat('adaptivity', 'prob_wx') if \
                config.has_option('adaptivity', 'prob_wx') else None, \
            opt_mthd_vx = config.get('adaptivity', 'opt_mthd_vx')[1:-1] if \
                config.has_option('adaptivity', 'opt_mthd_vx') else None, \
            opt_mthd_wx = config.get('adaptivity', 'opt_mthd_wx')[1:-1] if \
                config.has_option('adaptivity', 'opt_mthd_wx') else None, \
            blk_mks_r_v = eval(config.get('adaptivity', 'blk_mks_r_v')) if \
                config.has_option('adaptivity', 'blk_mks_r_v') else None, \
            blk_mks_r_w = eval(config.get('adaptivity', 'blk_mks_r_w')) if \
                config.has_option('adaptivity', 'blk_mks_r_w') else None
            )

    main_dict.update(aux_dict)

    return main_dict
    

def add_vx_blk_structure_info(parameters_str, n_blks, \
                                probblt, opt_mthd, blk_mks_r, blk_mks_v):

    """add information about the vx block structure to the parameters string"""

    config = ConfigParser.RawConfigParser(allow_no_value = True)
    config.readfp(BytesIO(parameters_str))

    no_adptv_section = False
    if not config.has_section('adaptivity'):
        config.add_section('adaptivity')
    else:
        config.remove_option('adaptivity', '/')

    blk_mks_r_str = ', '.join(map(str, blk_mks_r))
    blk_mks_v_str = ', '.join(map(str, blk_mks_v))

    # write adaptivity section
    config.set('adaptivity', 'n_vx_blks', n_blks)
    config.set('adaptivity', 'prob_vx', probblt)
    config.set('adaptivity', 'opt_mthd_vx',"'" + opt_mthd + "'")
    config.set('adaptivity', 'blk_mks_r_v', blk_mks_r_str)
    config.set('adaptivity', 'blk_mks_v', blk_mks_v_str)

    # modify parallel velocity range in the box section
    config.set('box', 'lv', blk_mks_v[0])
    
    configbytes = BytesIO(parameters_str)
    config.write(configbytes)
    exper_str = configbytes.getvalue()

    with BytesIO() as configbytes:
        config.write(configbytes)
        parameters_str = configbytes.getvalue()
        parameters_str += '\n/'

    # erase all empty line
    parameters_str = linesep.join([s for s in parameters_str.splitlines() if s])

    return parameters_str

def add_wx_blk_structure_info(parameters_str, n_blks, \
                                probblt, opt_mthd, blk_mks_r, blk_mks_w):

    """add information about the wx block structure to the parameters string"""

    config = ConfigParser.RawConfigParser(allow_no_value = True)
    config.readfp(BytesIO(parameters_str))

    no_adptv_section = False
    if not config.has_section('adaptivity'):
        config.add_section('adaptivity')
    else:
        config.remove_option('adaptivity', '/')

    blk_mks_r_str = ', '.join(map(str, blk_mks_r))
    blk_mks_w_str = ', '.join(map(str, blk_mks_w))

    # write adaptivity section
    config.set('adaptivity', 'n_wx_blks', n_blks)
    config.set('adaptivity', 'prob_wx', probblt)
    config.set('adaptivity', 'opt_mthd_wx',"'" + opt_mthd + "'")
    config.set('adaptivity', 'blk_mks_r_w', blk_mks_r_str)
    config.set('adaptivity', 'blk_mks_w', blk_mks_w_str)

    # modify parallel velocity range in the box section
    config.set('box', 'lw', blk_mks_w[0])
    
    configbytes = BytesIO(parameters_str)
    config.write(configbytes)
    exper_str = configbytes.getvalue()

    with BytesIO() as configbytes:
        config.write(configbytes)
        parameters_str = configbytes.getvalue()
        parameters_str += '\n/'

    # erase all empty line
    parameters_str = linesep.join([s for s in parameters_str.splitlines() if s])

    return parameters_str
