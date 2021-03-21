import sys
import json
import matplotlib.pyplot as plt
import numpy as np

times={
    "FE_Q":{
        "1proc":{
            "true":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            },
            "false":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            }
        },
        "2proc":{
            "true":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            },
            "false":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            }
        },
        "8proc":{
            "true":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            },
            "false":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            }
        },
        "16proc":{
            "true":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            },
            "false":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            }
        }
    },
    "FE_DGQ":{
        "1proc":{
            "true":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            },
            "false":{
                "combi10":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                },
                "combi100":{
                    "lmin1":[],
                    "lmin2":[],
                    "lmin3":[],
                    "fullgrids":[]
                }
            }
        }
    },
    "FE_Q_ngroup_4":{
        "lmin1":[],
        "lmin2":[],
        "lmin3":[],
        "fullgrids":[]
    }
}
dofs={
    "lmin1":[81,225,593,1505,3713,8961,21249,49665], #2 fehlen noch
    "lmin2":[125,425,1265,3489,9153,23169,57089,137729], #1 fehlt noch
    "lmin3":[729,2673,8289,23489,62849,161537,402945], #1 fehlt noch
    "fullgrids":[125,729,4913,35937,274625,2146689]
}
errors={
    "FE_Q":{
        "true":{
            "combi10":{
                "lmin1":[],
                "lmin2":[],
                "lmin3":[],
                "fullgrids":[]
            },
            "combi100":{
                "lmin1":[],
                "lmin2":[],
                "lmin3":[],
                "fullgrids":[]
            }
        },
        "false":{
            "combi10":{
                "lmin1":[],
                "lmin2":[],
                "lmin3":[],
                "fullgrids":[]
            },
            "combi100":{
                "lmin1":[],
                "lmin2":[],
                "lmin3":[],
                "fullgrids":[]
            }
        }
    },
    "FE_DGQ":{
        "true":{
            "combi10":{
                "lmin1":[],
                "lmin2":[],
                "lmin3":[],
                "fullgrids":[]
            },
            "combi100":{
                "lmin1":[],
                "lmin2":[],
                "lmin3":[],
                "fullgrids":[]
            }
        },
        "false":{
            "combi10":{
                "lmin1":[],
                "lmin2":[],
                "lmin3":[],
                "fullgrids":[]
            },
            "combi100":{
                "lmin1":[],
                "lmin2":[],
                "lmin3":[],
                "fullgrids":[]
            }
        }
    }
}
DG_dofs={
    "lmin1":[256,832,2432,6656,17408,44032,108544,262144],
    "lmin2":[512,2048,6656,19456,53248,139264,352256,868352],
    "lmin3":[4096,16384,53248,155648,425984,1114112,2818048],
    "fullgrids":[512,4096,32768,262144,2097152,16777216],
}
try:
    #loop over all attributes
    for method in ["Q","DGQ"]:
        for processes in [1,2,8,16]:
            if processes>1 and method=="DGQ":
                continue
            if processes==1:
                p="1,1,1"
            elif processes==2:
                p="1,2,1"
            elif processes==8:
                p="2,2,2"
            elif processes==16:
                p="4,2,2"

            for ncombi in [10,100]:
                dt=1.0/ncombi
                for DO_COMBINE in ["true","false"]:                
                    for lmax in range(2,10):   
                        s_lmax="_"+str(lmax)+","+str(lmax)+","+str(lmax)+"_"    
                        
                        for lmin in range(1,lmax+1):
                            if lmin>3:
                                if lmin<lmax:
                                    continue
                                elif lmax>7:
                                    continue
                            s_lmin="_"+str(lmin)+","+str(lmin)+","+str(lmin)+"_"
                            filepath='timers/helium/FE_'+method+DO_COMBINE+"/p_"+p+"_mi"+s_lmin+"ma"+s_lmax+"ngroup_1_ncombi_"+str(ncombi)+'.json'
                            try:
                                with open(filepath) as json_file:
                                    data = json.load(json_file)
                                    corr_comp=data['rank'+str(processes)]["events"]['correct computation'][0]
                                    lvl=""
                                    if lmin==lmax:
                                        lvl="fullgrids"
                                        if lmin<4:
                                            times["FE_"+method][str(processes)+"proc"][DO_COMBINE]["combi"+str(ncombi)]["lmin"+str(lmin)].append(corr_comp[1]-corr_comp[0])
                                    else:
                                        lvl="lmin"+str(lmin)
                                    times["FE_"+method][str(processes)+"proc"][DO_COMBINE]["combi"+str(ncombi)][lvl].append(corr_comp[1]-corr_comp[0])
                            except:
                                print(filepath)
                            

except:
    print(times)

for lmax in range(2,10):   
    s_lmax="_"+str(lmax)+","+str(lmax)+","+str(lmax)+"_"    
    
    for lmin in range(1,lmax+1):
        if lmin>3:
            if lmin<lmax:
                continue
            elif lmax>7:
                continue
        s_lmin="_"+str(lmin)+","+str(lmin)+","+str(lmin)+"_"
        filepath='timers/helium/FE_Qtrue/p_2,2,1_mi'+s_lmin+"ma"+s_lmax+"ngroup_4_ncombi_10.json"
        try:
            with open(filepath) as json_file:
                data = json.load(json_file)
                corr_comp=data['rank'+str(processes)]["events"]['correct computation'][0]
                lvl=""
                if lmin==lmax:
                    lvl="fullgrids"
                    if lmin<4:
                        times["FE_Q_ngroup_4"]["lmin"+str(lmin)].append(corr_comp[1]-corr_comp[0])
                else:
                    lvl="lmin"+str(lmin)
                times["FE_Q_ngroup_4"][lvl].append(corr_comp[1]-corr_comp[0])
        except:
            print(filepath)

#print(times)

f = open("L2_errors.dat", "r")
processes=1
for x in f:
    if not (x.startswith("Q ") or x.startswith("DGQ ")):
        name, val=x.split("=")
        method,do_combine,ncombi,min_lvl,max_lvl=name.split("_")
        lvl=""
        if processes==1:
            if min_lvl==max_lvl:
                lvl="fullgrids"
                if int(min_lvl)<4:
                    errors["FE_"+method][do_combine]["combi"+str(ncombi)]["lmin"+str(min_lvl)].append(float(val))
            else:
                lvl="lmin"+str(min_lvl)
            errors["FE_"+method][do_combine]["combi"+str(ncombi)][lvl].append(float(val))
    else:
        #compute the number of processes
        procs=x.split(", p=")[1]
        [one,two, three]=procs.split()
        print([one,two, three])
        if int(one)*int(two)*int(three)>1:
            processes=2
        else:
            processes=1
        
f.close()
#print(errors)

#plots error vs DOF decay
def CG_error_dof_10_true():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF")
    plt.plot(dofs["lmin1"],errors["FE_Q"]["true"]["combi10"]["lmin1"], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],errors["FE_Q"]["true"]["combi10"]["lmin2"], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(dofs["lmin3"],errors["FE_Q"]["true"]["combi10"]["lmin3"], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(dofs["fullgrids"],errors["FE_Q"]["true"]["combi10"]["fullgrids"], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error CG, ncombi=10, with combination")
    plt.legend()
    plt.savefig("tex_img/CG_error_dof_10_true.png")
    plt.close()

CG_error_dof_10_true()
#plots error vs DOF decay    
def CG_error_dof_100_true():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF")
    plt.plot(dofs["lmin1"],errors["FE_Q"]["true"]["combi100"]["lmin1"], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],errors["FE_Q"]["true"]["combi100"]["lmin2"], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(dofs["lmin3"],errors["FE_Q"]["true"]["combi100"]["lmin3"], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(dofs["fullgrids"],errors["FE_Q"]["true"]["combi100"]["fullgrids"], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error CG, ncombi=100, with combination")
    plt.legend()
    plt.savefig("tex_img/CG_error_dof_100_true.png")
    plt.close()

CG_error_dof_100_true()
#plots error vs DOF decay    
def CG_error_dof_10_false():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF")
    plt.plot(dofs["lmin1"],errors["FE_Q"]["false"]["combi10"]["lmin1"], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],errors["FE_Q"]["false"]["combi10"]["lmin2"], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(dofs["lmin3"],errors["FE_Q"]["false"]["combi10"]["lmin3"], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(dofs["fullgrids"],errors["FE_Q"]["false"]["combi10"]["fullgrids"], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error CG, ncombi=10, without combination")
    plt.legend()
    plt.savefig("tex_img/CG_error_dof_10_false.png")
    plt.close()

CG_error_dof_10_false()
#plots error vs DOF decay
def CG_error_dof_100_false():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF")
    plt.plot(dofs["lmin1"],errors["FE_Q"]["false"]["combi100"]["lmin1"], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],errors["FE_Q"]["false"]["combi100"]["lmin2"], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(dofs["lmin3"],errors["FE_Q"]["false"]["combi100"]["lmin3"], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(dofs["fullgrids"],errors["FE_Q"]["false"]["combi100"]["fullgrids"], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error CG, ncombi=100, without combination")
    plt.legend()
    plt.savefig("tex_img/CG_error_dof_100_false.png")
    plt.close()

CG_error_dof_100_false()
#plots error vs Time decay
def CG_error_times_10_true():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin1"],errors["FE_Q"]["true"]["combi10"]["lmin1"], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin2"],errors["FE_Q"]["true"]["combi10"]["lmin2"], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin3"],errors["FE_Q"]["true"]["combi10"]["lmin3"], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi10"]["fullgrids"],errors["FE_Q"]["true"]["combi10"]["fullgrids"], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error CG, ncombi=10, with combination")
    plt.legend()
    plt.savefig("tex_img/CG_error_times_10_true.png")
    plt.close()

CG_error_times_10_true()
#plots error vs time decay    
def CG_error_times_100_true():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin1"],errors["FE_Q"]["true"]["combi100"]["lmin1"], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin2"],errors["FE_Q"]["true"]["combi100"]["lmin2"], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin3"],errors["FE_Q"]["true"]["combi100"]["lmin3"], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi100"]["fullgrids"],errors["FE_Q"]["true"]["combi100"]["fullgrids"], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error CG, ncombi=100, with combination")
    plt.legend()
    plt.savefig("tex_img/CG_error_times_100_true.png")
    plt.close()

CG_error_times_100_true()
#plots error vs time decay
def CG_error_times_10_false():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi10"]["lmin1"],errors["FE_Q"]["false"]["combi10"]["lmin1"], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi10"]["lmin2"],errors["FE_Q"]["false"]["combi10"]["lmin2"], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi10"]["lmin3"],errors["FE_Q"]["false"]["combi10"]["lmin3"], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi10"]["fullgrids"],errors["FE_Q"]["false"]["combi10"]["fullgrids"], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error CG, ncombi=10, without combination")
    plt.legend()
    plt.savefig("tex_img/CG_error_times_10_false.png")
    plt.close()

CG_error_times_10_false()
#plots error vs time decay
def CG_error_times_100_false():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin1"],errors["FE_Q"]["false"]["combi100"]["lmin1"], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin2"],errors["FE_Q"]["false"]["combi100"]["lmin2"], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin3"],errors["FE_Q"]["false"]["combi100"]["lmin3"], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi100"]["fullgrids"],errors["FE_Q"]["false"]["combi100"]["fullgrids"], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error CG, ncombi=100, without combination")
    plt.legend()
    plt.savefig("tex_img/CG_error_times_100_false.png")
    plt.close()

CG_error_times_100_false()
#plots DG error vs time decay
def DG_error_times_10_false():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin1"][:6],errors["FE_DGQ"]["false"]["combi10"]["lmin1"][:6], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin2"][:6],errors["FE_DGQ"]["false"]["combi10"]["lmin2"][:6], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin3"][:5],errors["FE_DGQ"]["false"]["combi10"]["lmin3"][:5], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi10"]["fullgrids"][:5],errors["FE_DGQ"]["false"]["combi10"]["fullgrids"][:5], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error DG, ncombi=10, without combination")
    plt.legend()
    plt.savefig("tex_img/DG_error_times_10_false.png")
    plt.close()

DG_error_times_10_false()

def DG_error_dof_10_false():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF")
    plt.plot(DG_dofs["lmin1"][:6],errors["FE_DGQ"]["false"]["combi10"]["lmin1"][:6], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(DG_dofs["lmin2"][:6],errors["FE_DGQ"]["false"]["combi10"]["lmin2"][:6], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(DG_dofs["lmin3"][:5],errors["FE_DGQ"]["false"]["combi10"]["lmin3"][:5], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(DG_dofs["fullgrids"][:5],errors["FE_DGQ"]["false"]["combi10"]["fullgrids"][:5], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error DG, ncombi=10, without combination")
    plt.legend()
    plt.savefig("tex_img/DG_error_dofs_10_false.png")
    plt.close()

DG_error_dof_10_false()
#plots DG error vs time decay
def DG_error_times_10_true():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi10"]["lmin1"],errors["FE_DGQ"]["true"]["combi10"]["lmin1"][:1], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi10"]["lmin2"],errors["FE_DGQ"]["true"]["combi10"]["lmin2"][:2], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi10"]["lmin3"],errors["FE_DGQ"]["true"]["combi10"]["lmin3"][:1], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi10"]["fullgrids"][:5],errors["FE_DGQ"]["true"]["combi10"]["fullgrids"], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error DG, ncombi=10, with combination")
    plt.legend()
    plt.savefig("tex_img/DG_error_times_10_true.png")
    plt.close()

#DG_error_times_10_true()

def DG_error_dofs_10_true():
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(DG_dofs["lmin1"][:8],errors["FE_DGQ"]["true"]["combi10"]["lmin1"][:8], label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(DG_dofs["lmin2"][:7],errors["FE_DGQ"]["true"]["combi10"]["lmin2"][:7], label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(DG_dofs["lmin3"][:8],errors["FE_DGQ"]["true"]["combi10"]["lmin3"][:8], label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(DG_dofs["fullgrids"][:5],errors["FE_DGQ"]["true"]["combi10"]["fullgrids"][:5], label="Fullgrids",linestyle='dashdot', marker='v')
    plt.title("Error DG, ncombi=10, with combination")
    plt.legend()
    plt.savefig("tex_img/DG_error_dofs_10_true.png")
    plt.close()

DG_error_dofs_10_true()
#plots CG vs DG with combination and ncombi=10
def CG_vs_DG_10_true_times():
    print("ich vergeleiche CG und DG mit combi und 10 ncombi, times")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin1"],errors["FE_Q"]["true"]["combi10"]["lmin1"], label="CG, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin2"],errors["FE_Q"]["true"]["combi10"]["lmin2"], label="CG, $l_{min}=2$",linestyle='dashed', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin3"],errors["FE_Q"]["true"]["combi10"]["lmin3"], label="CG, $l_{min}=3$",linestyle='dotted', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi10"]["fullgrids"],errors["FE_Q"]["true"]["combi10"]["fullgrids"], label="CG, Fullgrids",linestyle='dashdot', marker='o')
    
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi10"]["lmin1"][:6],errors["FE_DGQ"]["true"]["combi10"]["lmin1"], label="DG, $l_{min}=1$",linestyle='solid', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi10"]["lmin2"][:6],errors["FE_DGQ"]["true"]["combi10"]["lmin2"], label="DG, $l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi10"]["lmin3"][:6],errors["FE_DGQ"]["true"]["combi10"]["lmin3"], label="DG, $l_{min}=3$",linestyle='dotted', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi10"]["fullgrids"][:5],errors["FE_DGQ"]["true"]["combi10"]["fullgrids"], label="DG, Fullgrids",linestyle='dashdot', marker='*')
    
    plt.title("CG vs DG, ncombi=10, with combination")
    plt.legend()
    plt.savefig("tex_img/CG_vs_DG_error_times_10_true.png")
    plt.close()

#CG_vs_DG_10_true_times()
#plots CG vs DG with combination and ncombi=10
def CG_vs_DG_10_true_dof():
    print("ich vergeleiche CG und DG mit combi und 10 ncombi, times")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF")
    plt.plot(dofs["lmin1"],errors["FE_Q"]["true"]["combi10"]["lmin1"], label="CG, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],errors["FE_Q"]["true"]["combi10"]["lmin2"], label="CG, $l_{min}=2$",linestyle='dashed', marker='o')
    plt.plot(dofs["lmin3"],errors["FE_Q"]["true"]["combi10"]["lmin3"], label="CG, $l_{min}=3$",linestyle='dotted', marker='o')
    plt.plot(dofs["fullgrids"],errors["FE_Q"]["true"]["combi10"]["fullgrids"], label="CG, Fullgrids",linestyle='dashdot', marker='o')
    
    plt.plot(dofs["lmin1"][:6],errors["FE_DGQ"]["true"]["combi10"]["lmin1"], label="DG, $l_{min}=1$",linestyle='solid', marker='*')
    plt.plot(dofs["lmin2"][:6],errors["FE_DGQ"]["true"]["combi10"]["lmin2"], label="DG, $l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(dofs["lmin3"][:5],errors["FE_DGQ"]["true"]["combi10"]["lmin3"], label="DG, $l_{min}=3$",linestyle='dotted', marker='*')
    plt.plot(dofs["fullgrids"][:5],errors["FE_DGQ"]["true"]["combi10"]["fullgrids"], label="DG, Fullgrids",linestyle='dashdot', marker='*')
    
    plt.title("CG vs DG, ncombi=10, with combination")
    plt.legend()
    plt.savefig("tex_img/CG_vs_DG_error_dofs_10_true.png")
    plt.close()

#CG_vs_DG_10_true_dof()
#plots CG vs DG with combination and ncombi=100
def CG_vs_DG_100_true_times():
    print("ich vergeleiche CG und DG mit combi und 100 ncombi, times")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin1"],errors["FE_Q"]["true"]["combi100"]["lmin1"], label="CG, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin2"],errors["FE_Q"]["true"]["combi100"]["lmin2"], label="CG, $l_{min}=2$",linestyle='dashed', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin3"],errors["FE_Q"]["true"]["combi100"]["lmin3"], label="CG, $l_{min}=3$",linestyle='dotted', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["true"]["combi100"]["fullgrids"],errors["FE_Q"]["true"]["combi100"]["fullgrids"], label="CG, Fullgrids",linestyle='dashdot', marker='o')
    
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi100"]["lmin1"],errors["FE_DGQ"]["true"]["combi100"]["lmin1"], label="DG, $l_{min}=1$",linestyle='solid', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi100"]["lmin2"],errors["FE_DGQ"]["true"]["combi100"]["lmin2"], label="DG, $l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi100"]["lmin3"],errors["FE_DGQ"]["true"]["combi100"]["lmin3"], label="DG, $l_{min}=3$",linestyle='dotted', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["true"]["combi100"]["fullgrids"],errors["FE_DGQ"]["true"]["combi100"]["fullgrids"], label="DG, Fullgrids",linestyle='dashdot', marker='*')
    
    plt.title("CG vs DG, ncombi=100, with combination")
    plt.legend()
    plt.savefig("tex_img/CG_vs_DG_error_times_100_true.png")
    plt.close()

#CG_vs_DG_100_true_times()
#plots CG vs DG with ncombi=10 and no combination
def CG_vs_DG_10_false_times():
    print("ich vergeleiche CG und DG mit combi und 10 ncombi, times")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi10"]["lmin1"],errors["FE_Q"]["false"]["combi10"]["lmin1"], label="CG, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi10"]["lmin2"],errors["FE_Q"]["false"]["combi10"]["lmin2"], label="CG, $l_{min}=2$",linestyle='dashed', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi10"]["lmin3"],errors["FE_Q"]["false"]["combi10"]["lmin3"], label="CG, $l_{min}=3$",linestyle='dotted', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi10"]["fullgrids"],errors["FE_Q"]["false"]["combi10"]["fullgrids"], label="CG, Fullgrids",linestyle='dashdot', marker='o')
    
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin1"][:6],errors["FE_DGQ"]["false"]["combi10"]["lmin1"][:6], label="DG, $l_{min}=1$",linestyle='solid', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin2"][:6],errors["FE_DGQ"]["false"]["combi10"]["lmin2"][:6], label="DG, $l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin3"][:5],errors["FE_DGQ"]["false"]["combi10"]["lmin3"][:5], label="DG, $l_{min}=3$",linestyle='dotted', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi10"]["fullgrids"][:5],errors["FE_DGQ"]["false"]["combi10"]["fullgrids"][:5], label="DG, Fullgrids",linestyle='dashdot', marker='*')
    
    plt.title("CG vs DG, ncombi=10, without combination")
    plt.legend()
    plt.savefig("tex_img/CG_vs_DG_error_times_10_false.png")
    plt.close()
#TODO:hier noch n besseren Vergleoich überlegen.
CG_vs_DG_10_false_times()

#plots CG vs DG with ncombi=10 and no combination
def CG_vs_DG_10_false_dof():
    print("ich vergeleiche CG und DG mit combi und 10 ncombi, times")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF'")
    plt.plot(dofs["lmin1"],errors["FE_Q"]["false"]["combi10"]["lmin1"], label="CG, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],errors["FE_Q"]["false"]["combi10"]["lmin2"], label="CG, $l_{min}=2$",linestyle='dashed', marker='o')
    plt.plot(dofs["lmin3"],errors["FE_Q"]["false"]["combi10"]["lmin3"], label="CG, $l_{min}=3$",linestyle='dotted', marker='o')
    plt.plot(dofs["fullgrids"],errors["FE_Q"]["false"]["combi10"]["fullgrids"], label="CG, Fullgrids",linestyle='dashdot', marker='o')
    
    plt.plot(DG_dofs["lmin1"][:6],errors["FE_DGQ"]["false"]["combi10"]["lmin1"][:6], label="DG, $l_{min}=1$",linestyle='solid', marker='*')
    plt.plot(DG_dofs["lmin2"][:6],errors["FE_DGQ"]["false"]["combi10"]["lmin2"][:6], label="DG, $l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(DG_dofs["lmin3"][:5],errors["FE_DGQ"]["false"]["combi10"]["lmin3"][:5], label="DG, $l_{min}=3$",linestyle='dotted', marker='*')
    plt.plot(DG_dofs["fullgrids"][:5],errors["FE_DGQ"]["false"]["combi10"]["fullgrids"][:5], label="DG, Fullgrids",linestyle='dashdot', marker='*')
    
    plt.title("CG vs DG, ncombi=10, without combination")
    plt.legend()
    plt.savefig("tex_img/CG_vs_DG_error_dof_10_false.png")
    plt.close()
#TODO: Diesen Vergleich dann auch hier einsetzen.
CG_vs_DG_10_false_dof()
#plots CG vs DG ncombi100 and no combination
def CG_vs_DG_100_false_times():
    print("ich vergeleiche CG und DG mit combi und 100 ncombi, times")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("Computation time")
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin1"],errors["FE_Q"]["false"]["combi100"]["lmin1"], label="CG, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin2"],errors["FE_Q"]["false"]["combi100"]["lmin2"], label="CG, $l_{min}=2$",linestyle='dashed', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin3"],errors["FE_Q"]["false"]["combi100"]["lmin3"], label="CG, $l_{min}=3$",linestyle='dotted', marker='o')
    plt.plot(times["FE_Q"]["1proc"]["false"]["combi100"]["fullgrids"],errors["FE_Q"]["false"]["combi100"]["fullgrids"], label="CG, Fullgrids",linestyle='dashdot', marker='o')
    
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi100"]["lmin1"],errors["FE_DGQ"]["false"]["combi100"]["lmin1"], label="DG, $l_{min}=1$",linestyle='solid', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi100"]["lmin2"],errors["FE_DGQ"]["false"]["combi100"]["lmin2"], label="DG, $l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi100"]["lmin3"],errors["FE_DGQ"]["false"]["combi100"]["lmin3"], label="DG, $l_{min}=3$",linestyle='dotted', marker='*')
    plt.plot(times["FE_DGQ"]["1proc"]["false"]["combi100"]["fullgrids"],errors["FE_DGQ"]["false"]["combi100"]["fullgrids"], label="DG, Fullgrids",linestyle='dashdot', marker='*')
    
    plt.title("CG vs DG, ncombi=100, without combination")
    plt.legend()
    plt.savefig("tex_img/CG_vs_DG_error_times_100_false.png")
    plt.close()

#CG_vs_DG_100_false_times()
#plots CG comparison of combination, ncombi=100
def CG_combi_compare_100():
    print("ich vergleiche CG_100_false mit CG_100_true")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF") 
    plt.plot(dofs["lmin1"],
    100*np.divide(np.subtract(errors["FE_Q"]["false"]["combi100"]["lmin1"],errors["FE_Q"]["true"]["combi100"]["lmin1"]),errors["FE_Q"]["false"]["combi100"]["lmin1"]), label="CG without combi, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],
    100*np.divide(np.subtract(errors["FE_Q"]["false"]["combi100"]["lmin2"],errors["FE_Q"]["true"]["combi100"]["lmin2"]),errors["FE_Q"]["false"]["combi100"]["lmin2"]), label="CG without combi, $l_{min}=2$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin3"],
    100*np.divide(np.subtract(errors["FE_Q"]["false"]["combi100"]["lmin3"],errors["FE_Q"]["true"]["combi100"]["lmin3"]),errors["FE_Q"]["false"]["combi100"]["lmin3"]), label="CG without combi, $l_{min}=3$",linestyle='solid', marker='o')
    
    plt.title("CG combination vs without combination, ncombi=100")
    plt.legend()
    plt.savefig("tex_img/CG_combi_compare_100_errors.png")
    plt.close()

    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("relative Time difference")
    plt.xlabel("#DOF") 
    plt.plot(dofs["lmin1"],
    100*np.divide(np.subtract(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin1"],times["FE_Q"]["1proc"]["false"]["combi100"]["lmin1"]),times["FE_Q"]["1proc"]["false"]["combi100"]["lmin1"]), label="CG without combi, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],
    100*np.divide(np.subtract(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin2"],times["FE_Q"]["1proc"]["false"]["combi100"]["lmin2"]),times["FE_Q"]["1proc"]["false"]["combi100"]["lmin2"]), label="CG without combi, $l_{min}=2$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin3"],
    100*np.divide(np.subtract(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin3"],times["FE_Q"]["1proc"]["false"]["combi100"]["lmin3"]),times["FE_Q"]["1proc"]["false"]["combi100"]["lmin3"]), label="CG without combi, $l_{min}=3$",linestyle='solid', marker='o')
    
    plt.title("CG combination vs without combination, ncombi=100")
    plt.legend()
    plt.savefig("tex_img/CG_combi_compare_100_times.png")
    plt.close()

CG_combi_compare_100()
#plots CG comparison of combination with ncombi=10
def CG_combi_compare_10():
    print("ich vergleiche CG_100_false mit CG_10_true")
    #plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF") 
    plt.ylim([0,0.5])
    plt.plot(dofs["lmin1"],
    np.divide(np.subtract(errors["FE_Q"]["false"]["combi10"]["lmin1"],errors["FE_Q"]["true"]["combi10"]["lmin1"]),errors["FE_Q"]["false"]["combi10"]["lmin1"]), label="$l_{min}=1$",linestyle='dashed', marker='o')
    plt.plot(dofs["lmin2"],
    np.divide(np.subtract(errors["FE_Q"]["false"]["combi10"]["lmin2"],errors["FE_Q"]["true"]["combi10"]["lmin2"]),errors["FE_Q"]["false"]["combi10"]["lmin2"]), label="$l_{min}=2$",linestyle='dotted', marker='*')
    plt.plot(dofs["lmin3"],
    np.divide(np.subtract(errors["FE_Q"]["false"]["combi10"]["lmin3"],errors["FE_Q"]["true"]["combi10"]["lmin3"]),errors["FE_Q"]["false"]["combi10"]["lmin3"]), label="$l_{min}=3$",linestyle='solid', marker='+')
    
    plt.title("CG combination vs without combination, ncombi=10")
    plt.legend()
    plt.savefig("tex_img/CG_combi_compare_10_errors.png")
    plt.close()

    #plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("relative Time difference")
    plt.xlabel("#DOF") 
    plt.ylim([-0.01,0.3])
    plt.plot(dofs["lmin1"],
    np.divide(np.subtract(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin1"],times["FE_Q"]["1proc"]["false"]["combi10"]["lmin1"]),times["FE_Q"]["1proc"]["false"]["combi10"]["lmin1"]), label="$l_{min}=1$",linestyle='dashed', marker='o')
    plt.plot(dofs["lmin2"],
    np.divide(np.subtract(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin2"],times["FE_Q"]["1proc"]["false"]["combi10"]["lmin2"]),times["FE_Q"]["1proc"]["false"]["combi10"]["lmin2"]), label="$l_{min}=2$",linestyle='dotted', marker='*')
    plt.plot(dofs["lmin3"],
    np.divide(np.subtract(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin3"],times["FE_Q"]["1proc"]["false"]["combi10"]["lmin3"]),times["FE_Q"]["1proc"]["false"]["combi10"]["lmin3"]), label="$l_{min}=3$",linestyle='solid', marker='+')
    
    plt.title("CG combination vs without combination, ncombi=10")
    plt.legend()
    plt.savefig("tex_img/CG_combi_compare_10_times.png")
    plt.close()

CG_combi_compare_10()
#plots DG comparison of combination with ncombi=100
def DG_combi_compare_100():
    print("ich vergleiche DG_100_false mit DG_10_true")
    #plt.yscale("log")
    plt.ylim([-1,1])
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF") 
    plt.plot(DG_dofs["lmin1"][:6],
    np.divide(np.subtract(errors["FE_DGQ"]["false"]["combi100"]["lmin1"][:6],errors["FE_DGQ"]["true"]["combi100"]["lmin1"][:6]),errors["FE_DGQ"]["false"]["combi100"]["lmin1"][:6]), label="DG without combi, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(DG_dofs["lmin2"][:6],
    np.divide(np.subtract(errors["FE_DGQ"]["false"]["combi100"]["lmin2"][:6],errors["FE_DGQ"]["true"]["combi100"]["lmin2"][:6]),errors["FE_DGQ"]["false"]["combi100"]["lmin2"][:6]), label="DG without combi, $l_{min}=2$",linestyle='solid', marker='o')
    plt.plot(DG_dofs["lmin3"][:5],
    np.divide(np.subtract(errors["FE_DGQ"]["false"]["combi100"]["lmin3"][:5],errors["FE_DGQ"]["true"]["combi100"]["lmin3"][:5]),errors["FE_DGQ"]["false"]["combi100"]["lmin3"][:5]), label="DG without combi, $l_{min}=3$",linestyle='solid', marker='o')
    
    plt.title("DG combination vs without combination, ncombi=100")
    plt.legend()
    plt.savefig("tex_img/DG_combi_compare_100_errors.png")
    plt.close()

    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("relative Time difference")
    plt.xlabel("#DOF") 
    # plt.plot(DG_dofs["lmin1"][:6],
    # np.divide(np.subtract(times["FE_DGQ"]["1proc"]["true"]["combi100"]["lmin1"][:6],times["FE_DGQ"]["1proc"]["false"]["combi100"]["lmin1"][:6]),times["FE_DGQ"]["1proc"]["false"]["combi100"]["lmin1"][:6]), label="DG without combi, $l_{min}=1$",linestyle='solid', marker='o')
    # plt.plot(DG_dofs["lmin2"][:6],
    # np.divide(np.subtract(times["FE_DGQ"]["1proc"]["true"]["combi100"]["lmin2"][:6],times["FE_DGQ"]["1proc"]["false"]["combi100"]["lmin2"][:6]),times["FE_DGQ"]["1proc"]["false"]["combi100"]["lmin2"][:6]), label="DG without combi, $l_{min}=2$",linestyle='solid', marker='o')
    # plt.plot(DG_dofs["lmin3"][:5],
    # np.divide(np.subtract(times["FE_DGQ"]["1proc"]["true"]["combi100"]["lmin3"][:5],times["FE_DGQ"]["1proc"]["false"]["combi100"]["lmin3"][:5]),times["FE_DGQ"]["1proc"]["false"]["combi100"]["lmin3"][:5]), label="DG without combi, $l_{min}=3$",linestyle='solid', marker='o')
    
    plt.title("DG combination vs without combination, ncombi=100")
    plt.legend()
    plt.savefig("tex_img/DG_combi_compare_100_times.png")
    plt.close()

#DG_combi_compare_100()
#plots DG comparison of combination with ncombi=10
def DG_combi_compare_10():
    print("ich vergleiche DG_10_false mit CG_10_true")
    #plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("Error")
    plt.xlabel("#DOF") 
    #plt.ylim([-100,1])
    index=4
    plt.plot(dofs["lmin1"][:index],
    np.divide(np.subtract(errors["FE_DGQ"]["false"]["combi10"]["lmin1"][:index],errors["FE_DGQ"]["true"]["combi10"]["lmin1"][:index]),errors["FE_DGQ"]["true"]["combi10"]["lmin1"][:index]), label="DG without combi, $l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"][:index],
    np.divide(np.subtract(errors["FE_DGQ"]["false"]["combi10"]["lmin2"][:index],errors["FE_DGQ"]["true"]["combi10"]["lmin2"][:index]),errors["FE_DGQ"]["true"]["combi10"]["lmin2"][:index]), label="DG without combi, $l_{min}=2$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin3"][:index],
    np.divide(np.subtract(errors["FE_DGQ"]["false"]["combi10"]["lmin3"][:index],errors["FE_DGQ"]["true"]["combi10"]["lmin3"][:index]),errors["FE_DGQ"]["true"]["combi10"]["lmin3"][:index]), label="DG without combi, $l_{min}=3$",linestyle='solid', marker='o')
    
    plt.title("DG combination vs without combination, ncombi=10")
    plt.legend()
    plt.savefig("tex_img/DG_combi_compare_10_errors.png")
    plt.close()

    #plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("relative Time difference")
    plt.xlabel("#DOF") 
    index=4
    plt.plot(dofs["lmin1"][:index],
    100*np.divide(np.subtract(times["FE_DGQ"]["1proc"]["true"]["combi10"]["lmin1"][:index],times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin1"][:index]),times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin1"][:index]), label="DG without combi, $l_{min}=1$",linestyle='solid', marker='o')
    index=2
    plt.plot(dofs["lmin2"][:index],
    100*np.divide(np.subtract(times["FE_DGQ"]["1proc"]["true"]["combi10"]["lmin2"][:index],times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin2"][:index]),times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin2"][:index]), label="DG without combi, $l_{min}=2$",linestyle='solid', marker='o')
    index=4
    plt.plot(dofs["lmin3"][:index],
    100*np.divide(np.subtract(times["FE_DGQ"]["1proc"]["true"]["combi10"]["lmin3"][:index],times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin3"][:index]),times["FE_DGQ"]["1proc"]["false"]["combi10"]["lmin3"][:index]), label="DG without combi, $l_{min}=3$",linestyle='solid', marker='o')
    
    plt.title("DG combination vs without combination, ncombi=10")
    plt.legend()
    plt.savefig("tex_img/DG_combi_compare_10_times.png")
    plt.close()
    

DG_combi_compare_10()
def parallelization_speedup():
    print("Ich vergleiche ")
    plt.xscale("log")
    plt.ylabel("Computation time")
    plt.yscale("log")
    plt.xlabel("#processes")
    plt.plot([1,2,8],[times["FE_Q"]["1proc"]["true"]["combi10"]["lmin1"][7],times["FE_Q"]["2proc"]["true"]["combi10"]["lmin1"][7],
    times["FE_Q"]["8proc"]["true"]["combi10"]["lmin1"][7]],label="$l_{\mathrm{min}}=1, l_{\mathrm{max}}=9$",marker='o', linestyle='dashed')
    
    plt.plot([1,2,8,16],[times["FE_Q"]["1proc"]["true"]["combi10"]["lmin2"][6],times["FE_Q"]["2proc"]["true"]["combi10"]["lmin2"][6],
    times["FE_Q"]["8proc"]["true"]["combi10"]["lmin2"][6],times["FE_Q"]["16proc"]["true"]["combi10"]["lmin2"][6]],label="$l_{\mathrm{min}}=2, l_{\mathrm{max}}=8$",linestyle='solid',marker='+')
    
    plt.plot([1,2,8,16],[times["FE_Q"]["1proc"]["true"]["combi10"]["fullgrids"][5],times["FE_Q"]["2proc"]["true"]["combi10"]["fullgrids"][5],
    times["FE_Q"]["8proc"]["true"]["combi10"]["fullgrids"][5],times["FE_Q"]["16proc"]["true"]["combi10"]["fullgrids"][5]],label="Fullgrid l=7",linestyle='dotted',marker='*')
    plt.plot([16], [times["FE_Q_ngroup_4"]["lmin1"][4]],label="4 Processing groups, $l_{\mathrm{min}}=1, l_{\mathrm{max}}=9$",linestyle='solid',marker='+')
    plt.plot([16], [times["FE_Q_ngroup_4"]["lmin2"][2]],label="4 Processing groups, $l_{\mathrm{min}}=2, l_{\mathrm{max}}=8$",linestyle='solid',marker='+')
    plt.legend()
    plt.xticks([1,2,8,16],[1,2,8,16])
    plt.savefig("tex_img/parallel.png")
    plt.close()



    # plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["false"]["combi10"]["lmin1"],times["FE_Q"]["2proc"]["false"]["combi10"]["lmin1"]), label="2 proc",linestyle='solid', marker='o')
    # plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["false"]["combi10"]["lmin1"],times["FE_Q"]["8proc"]["false"]["combi10"]["lmin1"]), label="8 proc",linestyle='solid', marker='o')
    # #plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["false"]["combi10"]["lmin1"],times["FE_Q"]["16proc"]["false"]["combi10"]["lmin1"]), label="16 proc",linestyle='solid', marker='o')
    # plt.title("Speedup by using multiple processes for lmin=1, dt=0.1, no combi")
    # plt.legend()
    # plt.savefig("tex_img/parallel.png")
    # plt.close()

    # plt.xscale("log")
    # plt.ylabel("Speedup")
    # plt.xlabel("#DOF")
    # plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin1"],times["FE_Q"]["2proc"]["false"]["combi100"]["lmin1"]), label="2 proc",linestyle='solid', marker='o')
    # plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin1"],times["FE_Q"]["8proc"]["false"]["combi100"]["lmin1"]), label="8 proc",linestyle='solid', marker='o')
    # #plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin1"],times["FE_Q"]["16proc"]["false"]["combi100"]["lmin1"]), label="16 proc",linestyle='solid', marker='o')
    # plt.title("Speedup by using multiple processes for lmin=1, dt=0.01, no combi")
    # plt.legend()
    # plt.savefig("tex_img/parallel2.png")
    # plt.close()

    # plt.xscale("log")
    # plt.ylabel("Speedup")
    # plt.xlabel("#DOF")
    # plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin1"],times["FE_Q"]["2proc"]["true"]["combi100"]["lmin1"]), label="2 proc",linestyle='solid', marker='o')
    # plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin1"],times["FE_Q"]["8proc"]["true"]["combi100"]["lmin1"]), label="8 proc",linestyle='solid', marker='o')
    # #plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["true"]["combi100"]["lmin1"],times["FE_Q"]["16proc"]["true"]["combi100"]["lmin1"]), label="16 proc",linestyle='solid', marker='o')
    # plt.title("Speedup by using multiple processes for lmin=1, dt=0.01, with combi")
    # plt.legend()
    # plt.savefig("tex_img/parallel3.png")
    # plt.close()

    # plt.xscale("log")
    # plt.ylabel("Speedup")
    # plt.xlabel("#DOF")
    # plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin1"],times["FE_Q"]["2proc"]["true"]["combi10"]["lmin1"]), label="2 proc",linestyle='solid', marker='o')
    # plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin1"],times["FE_Q"]["8proc"]["true"]["combi10"]["lmin1"]), label="8 proc",linestyle='solid', marker='o')
    # #plt.plot(dofs["lmin1"],np.divide(times["FE_Q"]["1proc"]["true"]["combi10"]["lmin1"],times["FE_Q"]["16proc"]["true"]["combi10"]["lmin1"]), label="16 proc",linestyle='solid', marker='o')
    # plt.title("Speedup by using multiple processes for lmin=1, dt=0.1, with combi")
    # plt.legend()
    # plt.savefig("tex_img/parallel4.png")
    # plt.close()

parallelization_speedup()
#compares ncombi10 to ncombi100
def CG_dt_compare_false():
    print("Ich bin ein Vergleich zwischen den Zeitschrittweiten. Bringt es mehr zeitschritte zu berechnen oder reichen 10 aus?")
    #plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("relative error")
    plt.xlabel("#DOF")
    plt.ylim([0,0.042])
    plt.yticks([0,0.02,0.04],[0,0.02,0.04])
    plt.plot(dofs["lmin1"],np.divide(np.subtract(errors["FE_Q"]["true"]["combi10"]["lmin1"],errors["FE_Q"]["true"]["combi100"]["lmin1"]),errors["FE_Q"]["true"]["combi10"]["lmin1"]), label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],np.divide(np.subtract(errors["FE_Q"]["true"]["combi10"]["lmin2"],errors["FE_Q"]["true"]["combi100"]["lmin2"]),errors["FE_Q"]["true"]["combi10"]["lmin2"]), label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(dofs["lmin3"][:6],np.divide(np.subtract(errors["FE_Q"]["true"]["combi10"]["lmin3"][:6],errors["FE_Q"]["true"]["combi100"]["lmin3"][:6]),errors["FE_Q"]["true"]["combi10"]["lmin3"][:6]), label="$l_{min}=3$",linestyle='dotted', marker='+')
    #plt.plot(dofs["fullgrids"],np.abs(np.subtract(errors["FE_Q"]["true"]["combi10"]["fullgrids"],errors["FE_Q"]["true"]["combi100"]["fullgrids"])), label="Fullgrids",linestyle='dashdot', marker='v')
    #plt.title("Absolute difference of CG ncombi=10 to CG with ncombi=100, with combination.")
    plt.legend()
    plt.savefig("tex_img/CG_error_dof_compare_false.png")
    plt.close()
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("relative time difference")
    plt.xlabel("#DOF")
    plt.yticks([0.3,1,2,3,4],[0.3,1,2,3,4])
    plt.plot(dofs["lmin1"],np.divide(np.subtract(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin1"],times["FE_Q"]["1proc"]["false"]["combi10"]["lmin1"]),times["FE_Q"]["1proc"]["false"]["combi10"]["lmin1"]), label="$l_{min}=1$",linestyle='solid', marker='o')
    plt.plot(dofs["lmin2"],np.divide(np.subtract(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin2"],times["FE_Q"]["1proc"]["false"]["combi10"]["lmin2"]),times["FE_Q"]["1proc"]["false"]["combi10"]["lmin2"]), label="$l_{min}=2$",linestyle='dashed', marker='*')
    plt.plot(dofs["lmin3"],np.divide(np.subtract(times["FE_Q"]["1proc"]["false"]["combi100"]["lmin3"],times["FE_Q"]["1proc"]["false"]["combi10"]["lmin3"]),times["FE_Q"]["1proc"]["false"]["combi10"]["lmin3"]), label="$l_{min}=3$",linestyle='dotted', marker='+')
    plt.plot(dofs["fullgrids"],np.divide(np.subtract(times["FE_Q"]["1proc"]["false"]["combi100"]["fullgrids"],times["FE_Q"]["1proc"]["false"]["combi10"]["fullgrids"]),times["FE_Q"]["1proc"]["false"]["combi10"]["fullgrids"]), label="Fullgrids",linestyle='dashdot', marker='v')
    #plt.title("Absolute difference of CG ncombi=10 to CG with ncombi=100, without combination.")
    plt.legend()
    plt.savefig("tex_img/CG_times_dof_compare_false.png")
    plt.close()

CG_dt_compare_false()
