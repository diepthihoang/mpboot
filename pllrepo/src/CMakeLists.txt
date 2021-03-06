#set( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -msse3 -DRAXML_USE_LLALLOC -D_USE_PTHREADS -D_OPTIMIZED_FUNCTIONS -D__SIM_SSE3 -fno-builtin" )

#add_executable( raxml_light axml.c  optimizeModel.c trash.c searchAlgo.c topologies.c treeIO.c models.c evaluatePartialGenericSpecial.c evaluateGenericSpecial.c newviewGenericSpecial.c makenewzGenericSpecial.c bipartitionList.c restartHashTable.c fastDNAparsimony.c randomTree.c lockless_allocator/ll_alloc.c mem_alloc.c recom.c)

#target_link_libraries( raxml_light m pthread )

if ((IQTREE_FLAGS MATCHES "avx" OR IQTREE_FLAGS MATCHES "fma") AND IQTREE_FLAGS MATCHES "omp")
  add_library(pll 
  avxLikelihood.c
  alignment.c
  bipartitionList.c
  evaluateGenericSpecial.c
  evaluatePartialGenericSpecial.c
  fastDNAparsimony.c
  hardware.c
  hash.c
  lexer.c
  makenewzGenericSpecial.c
  models.c
  newick.c
  newviewGenericSpecial.c
  genericParallelization.c
  optimizeModel.c
  parsePartition.c
  queue.c
  randomTree.c
  recom.c
  restartHashTable.c
  searchAlgo.c
  ssort.c
  stack.c
  topologies.c
  trash.c
  treeIO.c
  utils.c)
elseif (IQTREE_FLAGS MATCHES "omp")
  add_library(pll 
  alignment.c
  bipartitionList.c
  evaluateGenericSpecial.c
  evaluatePartialGenericSpecial.c
  fastDNAparsimony.c
  hardware.c
  hash.c
  lexer.c
  makenewzGenericSpecial.c
  models.c
  newick.c
  genericParallelization.c
  newviewGenericSpecial.c
  optimizeModel.c
  parsePartition.c
  queue.c
  randomTree.c
  recom.c
  restartHashTable.c
  searchAlgo.c
  ssort.c
  stack.c
  topologies.c
  trash.c
  treeIO.c
  utils.c)
elseif (IQTREE_FLAGS MATCHES "avx" OR IQTREE_FLAGS MATCHES "fma")
  add_library(pll 
    avxLikelihood.c
    alignment.c
    bipartitionList.c
    evaluateGenericSpecial.c
    evaluatePartialGenericSpecial.c
    fastDNAparsimony.c
    hardware.c
    hash.c
    lexer.c
    makenewzGenericSpecial.c
    models.c
    newick.c
    newviewGenericSpecial.c
    optimizeModel.c
    parsePartition.c
    queue.c
    randomTree.c
    recom.c
    restartHashTable.c
    searchAlgo.c
    ssort.c
    stack.c
    topologies.c
    trash.c
    treeIO.c
    utils.c)
else ()
  add_library(pll 
    alignment.c
    bipartitionList.c
    evaluateGenericSpecial.c
    evaluatePartialGenericSpecial.c
    fastDNAparsimony.c
    hardware.c
    hash.c
    lexer.c
    makenewzGenericSpecial.c
    models.c
    newick.c
    newviewGenericSpecial.c
    optimizeModel.c
    parsePartition.c
    queue.c
    randomTree.c
    recom.c
    restartHashTable.c
    searchAlgo.c
    ssort.c
    stack.c
    topologies.c
    trash.c
    treeIO.c
    utils.c)       
endif()
