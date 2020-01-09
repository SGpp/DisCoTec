#!/usr/bin/env python
""" Python script for plotting of profile data

 Using the classes ProfileData and PlotNeodata,
 as well as FSAmomData
"""

import argparse

import ProfileData
import FSAmomData
import PlotNeodata as Pnd
# This set controls what is plotted.
# Without the aprropriate element here, the plotting routines return nothing.
PLOTTHIS = [
    'plotT',
    'plotchi',
    'plotQ',
    'plotch',
    'plotlocQ',
#    'plotP',
    'plotj',
#    'plotG',
#    'plotturb',
    'plotnu',
    'plotfsa',
    'ploterad'
    ]
PLOTTHIS = set(PLOTTHIS)


parser = argparse.ArgumentParser()
parser.add_argument("time", help="Time window for averaging in the form start:end. "
                    "Use f for the first and l for the last time in the profile file.")
parser.add_argument("fileextension", nargs='+', help="List of file extensions to be processed.")
parser.add_argument("--turb", "-t", action='store_true',
                    help="Include turbulent transport in plots")
parser.add_argument("--xtcont", "-x", action='store_true',
                    help="Create x-t contour plots instead of time averaged profiles")
parser.add_argument("--contin", "-c", action='store_true',
                    help="Try to combine the given files into one averaged profile"
                    "This requires matching parameter sets")
parser.add_argument("--makepdf", "-p", action='store_true', help="Generate pdf files of the plots")

args = parser.parse_args()

fileextension = args.fileextension

if args.turb:
    PLOTTHIS.add('plotturb')

if args.time.split(':')[0] == 'f':
    starttime = -1
elif args.time.split(':')[0] == 'l':
    starttime = -2
else:
    starttime = int(float(args.time.split(':')[0]))
if args.time.split(':')[1] == 'f':
    endtime = -1
elif args.time.split(':')[1] == 'l':
    endtime = -2
else:
    endtime = int(float(args.time.split(':')[1]))

inputdata = []
inputfsa = []

for fext in fileextension:
    prdat = ProfileData.ProfileData(fext, starttime, endtime, PLOTTHIS)
    fsadat = FSAmomData.FSAmomData(fext, starttime, endtime, PLOTTHIS)
    prdat.get_profiles()
    fsadat.getfsadata()
    inputdata.append(prdat)
    inputfsa.append(fsadat)
if args.contin:
        inputdata = [ProfileData.glueprofiledata(inputdata)]
#            del(inputdata[i])

# inputdata[1].plotset.remove('plotch')
if args.xtcont:
    plot = Pnd.PlotNC_xt(inputdata)
    plot.plotT_cont(args.makepdf)
    plot.plotQ_cont(args.makepdf)
    plot.plotjbs_cont(args.makepdf)
    plot.show()
else:
    plot = Pnd.PlotNeodata(inputdata, inputfsa)
    plot.plot_T(args.makepdf)
    plot.plot_Q(args.makepdf)
    plot.plot_chi(args.makepdf)
    plot.plot_jbs(args.makepdf)
    plot.plot_nustar(args.makepdf)
    plot.plot_erad(args.makepdf)
    plot.plot_fsa(args.makepdf)
    plot.show()
