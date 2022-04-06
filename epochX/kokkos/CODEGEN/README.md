# madgraph_kokkos_plugin
A Plugin for Kokkos in MG

# install instructions

```bash

# checkout madgraph GPU version
bzr branch lp:~maddevelopers/mg5amcnlo/2.7.0_gpu

# copy all the codegen files to madgraph:
cd whereYouCheckedOut_MG2.7.0
rm -rf PLUGIN/*; cp -dpr ../madevent_gpu/epochX/kokkos/CODEGEN/PLUGIN/KOKKOS_SA_OUTPUT ./PLUGIN/

# make sure you delete the old generated code if it exists
rm -rf CODEGEN_kokkos_gg_tt

# Generate the code
echo -e "generate g g > t t~\noutput standalone_kokkos CODEGEN_kokkos_gg_tt\n"|./bin/mg5_aMC

```
