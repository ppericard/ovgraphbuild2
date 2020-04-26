cmake_minimum_required (VERSION 3.4)
project (ovgraphbuild2 CXX)

find_package (SeqAn3 3.0.0 REQUIRED HINTS "${CMAKE_SOURCE_DIR}/../seqan3/build_system")

add_executable (ovgraphbuild ovgraphbuild.cpp)

target_link_libraries (ovgraphbuild seqan3::seqan3)
