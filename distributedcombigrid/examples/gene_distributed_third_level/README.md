Running the third-level gene example:
-------------------------------------

Under `distributedcombigrid/examples/gene_distributed_third_level`  you can find
a third-level combination example with a Gene task. The following instructions
explain how to run the example.

Prerequisites are that the compilation of distributedcombigrid, the third-level
manager and the example (requires CombiDataType to be complex which can be set
in `distributedcombigrid/src/sgpp/distributedcombigrid/utils/Config.hpp`) was 
successful.

1. The setup of the example depends on the `gene_python_interface`. If there is
   no corresponding folder under `combi/` ask someone who has access and copy it
   there.

2. Create a parameter file that suits your setup (Use ctparam_tl as a template)
   and set the variable `paramfile` accordingly in `preproc.py`.

3. Run the `third-level-manager` which must be accessible form both systems.

4. Start `run.sh` on both systems.


If you plan to run the example on hazelhen there are a few more things to
mention:

1. Due to firewall restrictions, all connections must be tunneled with ssh and
   only outgoing connections are allowed on the compute nodes. 

2. Let's suppose the third level manager and rabbitmq server listen on port 9999
   at neon.

Use the following command on neon to create a tunnel which allows accessing both
servers from the login node:

`ssh -R 9999:localhost:9999 xy@eslogin002.hww.de `

3. Sadly, the ssh server configuration on the login nodes does not allow direct
   access from the compute nodes to the ssh server now sitting on port 9999 at
   the login node. To circumvent this restriction we can set up a temporary port
   forwarding from another port to 9999. Therefore you can use the `forward.sh`
   script located under ... . 
