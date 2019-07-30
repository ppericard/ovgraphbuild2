cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)

set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package (SeqAn3 REQUIRED)

add_executable (ovgraphbuild ovgraphbuild.cpp)

target_link_libraries (ovgraphbuild seqan3::seqan3)
