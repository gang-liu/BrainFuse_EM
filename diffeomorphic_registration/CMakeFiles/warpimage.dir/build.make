# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/cmake/2.8.9/bin/cmake

# The command to remove a file.
RM = /autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/cmake/2.8.9/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/cmake/2.8.9/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration

# Include any dependencies generated for this target.
include CMakeFiles/warpimage.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/warpimage.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/warpimage.dir/flags.make

CMakeFiles/warpimage.dir/mex_warpimage.cpp.o: CMakeFiles/warpimage.dir/flags.make
CMakeFiles/warpimage.dir/mex_warpimage.cpp.o: mex_warpimage.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/warpimage.dir/mex_warpimage.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -D_GNU_SOURCE -pthread -D_FILE_OFFSET_BITS=64 -DMX_COMPAT_32 -o CMakeFiles/warpimage.dir/mex_warpimage.cpp.o -c /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration/mex_warpimage.cpp

CMakeFiles/warpimage.dir/mex_warpimage.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/warpimage.dir/mex_warpimage.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -D_GNU_SOURCE -pthread -D_FILE_OFFSET_BITS=64 -DMX_COMPAT_32 -E /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration/mex_warpimage.cpp > CMakeFiles/warpimage.dir/mex_warpimage.cpp.i

CMakeFiles/warpimage.dir/mex_warpimage.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/warpimage.dir/mex_warpimage.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -fPIC -D_GNU_SOURCE -pthread -D_FILE_OFFSET_BITS=64 -DMX_COMPAT_32 -S /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration/mex_warpimage.cpp -o CMakeFiles/warpimage.dir/mex_warpimage.cpp.s

CMakeFiles/warpimage.dir/mex_warpimage.cpp.o.requires:
.PHONY : CMakeFiles/warpimage.dir/mex_warpimage.cpp.o.requires

CMakeFiles/warpimage.dir/mex_warpimage.cpp.o.provides: CMakeFiles/warpimage.dir/mex_warpimage.cpp.o.requires
	$(MAKE) -f CMakeFiles/warpimage.dir/build.make CMakeFiles/warpimage.dir/mex_warpimage.cpp.o.provides.build
.PHONY : CMakeFiles/warpimage.dir/mex_warpimage.cpp.o.provides

CMakeFiles/warpimage.dir/mex_warpimage.cpp.o.provides.build: CMakeFiles/warpimage.dir/mex_warpimage.cpp.o

CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o: CMakeFiles/warpimage.dir/flags.make
CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o: /autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c
	$(CMAKE_COMMAND) -E cmake_progress_report /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -fPIC -D_GNU_SOURCE -pthread -D_FILE_OFFSET_BITS=64 -DMX_COMPAT_32 -o CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o   -c /autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c

CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -fPIC -D_GNU_SOURCE -pthread -D_FILE_OFFSET_BITS=64 -DMX_COMPAT_32 -E /autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c > CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.i

CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -fPIC -D_GNU_SOURCE -pthread -D_FILE_OFFSET_BITS=64 -DMX_COMPAT_32 -S /autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c -o CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.s

CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o.requires:
.PHONY : CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o.requires

CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o.provides: CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o.requires
	$(MAKE) -f CMakeFiles/warpimage.dir/build.make CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o.provides.build
.PHONY : CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o.provides

CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o.provides.build: CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o

# Object files for target warpimage
warpimage_OBJECTS = \
"CMakeFiles/warpimage.dir/mex_warpimage.cpp.o" \
"CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o"

# External object files for target warpimage
warpimage_EXTERNAL_OBJECTS =

warpimage.mexa64: CMakeFiles/warpimage.dir/mex_warpimage.cpp.o
warpimage.mexa64: CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o
warpimage.mexa64: CMakeFiles/warpimage.dir/build.make
warpimage.mexa64: /autofs/space/lyon_006/pubsw/common/matlab/7.4/bin/glnxa64/libmx.so
warpimage.mexa64: /autofs/space/lyon_006/pubsw/common/matlab/7.4/bin/glnxa64/libmex.so
warpimage.mexa64: /autofs/space/lyon_006/pubsw/common/matlab/7.4/bin/glnxa64/libmat.so
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libITKCommon-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libITKStatistics-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libITKCommon-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libitksys-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libITKVNLInstantiation-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libitkvnl_algo-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libitkv3p_lsqr-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libitkvnl-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libitkvcl-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libitkdouble-conversion-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libitkNetlibSlatec-4.5.so.1
warpimage.mexa64: /autofs/cluster/con_003/users/mert/lukeliu/Download/debug_ITKandStatismo/build_itkforBFL/lib/libitkv3p_netlib-4.5.so.1
warpimage.mexa64: CMakeFiles/warpimage.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library warpimage.mexa64"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/warpimage.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/warpimage.dir/build: warpimage.mexa64
.PHONY : CMakeFiles/warpimage.dir/build

CMakeFiles/warpimage.dir/requires: CMakeFiles/warpimage.dir/mex_warpimage.cpp.o.requires
CMakeFiles/warpimage.dir/requires: CMakeFiles/warpimage.dir/autofs/space/lyon_006/pubsw/common/matlab/7.4/extern/src/mexversion.c.o.requires
.PHONY : CMakeFiles/warpimage.dir/requires

CMakeFiles/warpimage.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/warpimage.dir/cmake_clean.cmake
.PHONY : CMakeFiles/warpimage.dir/clean

CMakeFiles/warpimage.dir/depend:
	cd /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration /autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_test/diffeomorphic_registration/CMakeFiles/warpimage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/warpimage.dir/depend
