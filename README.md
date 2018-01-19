# mpboot
MPBoot: Fast phylogenetic maximum parsimony tree inference and bootstrap approximation

COMPILING INSTRUCTION
* Clone the source code, unzip it, and rename to "source"
* Create folder "build" outside folder "source"
* Change directory to "build"
* Run cmake with command: cmake ../source -DIQTREE_FLAGS=sse4
* Replace "sse4" by "avx" in above command if you decide to run MPBoot on AVX architecture
* Run make
* You will find the executable name "mpboot" once the make command is done.
