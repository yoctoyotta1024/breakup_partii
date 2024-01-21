To build, run and plot entire experiments:

note:
coalbu = coalescence + breakup
coalre = coalescence + rebound
coalbure = coalescence + breakup + rebound
coalonly = coalescence only (original Long kernel)

1) choose which cc[X]_config_collisions.hpp to copy into config_collision.hpp. E.g. for experiment cc1:
cd src/src/
cp cc1_config_collisions.hpp config_collisions.hpp

2) build experiments:
cd src/
./build_buii.sh

3) make initial conditions using wet or dry python scripts. E.g. for wet script:
cd src/
./initSDs.sh wet 

4) execute runs for each experiment:
cd src/
./run_buii.sh

5) make ensemble dataset from runs for each experiment:
cd src_ensemb/
./ensemb_buii.sh

6) optional plotting of data for a single run of an experiment:
cd quickplt/
./quickplotrun_buii.sh [experiment_name] [run_number]
e.g. ./quickplotrun_buii.sh coalbure 0

7) optional plotting of ensemble data for an experiment:
cd quickplt/
./quickplotens_buii.sh [experiment_name]

8) optional plotting of ensemble data from several experiments:
cd ensplt/
./pltens_buii.sh

