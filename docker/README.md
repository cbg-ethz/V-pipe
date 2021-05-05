# Docker setup for vpipe

1. Docker images are named vpipe-${BRANCHNAME}

2. This folder holds scripts which can run from any working directory.

3. To bootstrap a project:

   $ mkdir work
   $ cd work
   $ ${PATHTOTHISFOLDER}/init_project.sh

   this will create vpipe.config and references/

4. To run vpipe:
   $ cd work
   $ mkdir samples
   $ # populate samples folder
   $ # maybe edit vpipe-config
   $ ${PATHTOTHISFOLDER}/vpipe_docker.sh -p -j 2

   Technical note: the .snakemake folder is mounted into the current working directory
   on the host. Thus conda packages and conda environments are not created within
   the container but on the host file system.
   This reduces container size and persists state between different calls of 
   vpipe_docker.sh
