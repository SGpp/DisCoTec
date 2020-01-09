#include <mpi.h>
#include <stdio.h>
#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>

int main(int argc, char *argv[]) {

  _BGP_Personality_t personality;
  int xnodes, ynodes, znodes, coreID, myX,myY,myZ;
  int pset_size, node_config, taskid, ntasks, itask;
  char location[128];

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);

  Kernel_GetPersonality(&personality, sizeof(personality));

  if ( taskid == 0 ) {
    xnodes = personality.Network_Config.Xnodes;
    ynodes = personality.Network_Config.Ynodes;
    znodes = personality.Network_Config.Znodes;
    node_config = personality.Kernel_Config.ProcessConfig;
    pset_size = personality.Network_Config.PSetSize;
    /*pset_rank = personality.Network_Config.RankInPSet;*/

    printf("torus dimensions = <%d,%d,%d>\n",xnodes, ynodes,znodes);
    if (node_config == _BGP_PERS_PROCESSCONFIG_SMP) printf("SMP mode\n");
    else if (node_config == _BGP_PERS_PROCESSCONFIG_VNM) printf("VN mode\n");
    else if (node_config == _BGP_PERS_PROCESSCONFIG_2x2) printf("DUAL mode\n");
    else printf("Unkown mode.\n");

    printf("Number of processors in pset = %d, no. of MPI tasks = %d\n",pset_size,ntasks);

  }

  MPI_Barrier(MPI_COMM_WORLD);
  myX = personality.Network_Config.Xcoord;
  myY = personality.Network_Config.Ycoord;
  myZ = personality.Network_Config.Zcoord;

  coreID = Kernel_PhysicalProcessorID();
  BGP_Personality_getLocationString(&personality, location);
  for (itask=0;itask<ntasks;itask++) {
    if ( itask == taskid ) printf("MPI rank %d has torus coords <%d,%d,%d>  cpu = %d,        location = %s\n", taskid, myX, myY, myZ, coreID, location);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Finalize();
}
