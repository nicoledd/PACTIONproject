# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.27.4/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.27.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/nsdong2/Downloads/lemon-1.3.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/nsdong2/Downloads/lemon-1.3.1/build

# Include any dependencies generated for this target.
include tools/CMakeFiles/dimacs-to-lgf.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tools/CMakeFiles/dimacs-to-lgf.dir/compiler_depend.make

# Include the progress variables for this target.
include tools/CMakeFiles/dimacs-to-lgf.dir/progress.make

# Include the compile flags for this target's objects.
include tools/CMakeFiles/dimacs-to-lgf.dir/flags.make

tools/CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.o: tools/CMakeFiles/dimacs-to-lgf.dir/flags.make
tools/CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.o: /Users/nsdong2/Downloads/lemon-1.3.1/tools/dimacs-to-lgf.cc
tools/CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.o: tools/CMakeFiles/dimacs-to-lgf.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/nsdong2/Downloads/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tools/CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.o"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/tools && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tools/CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.o -MF CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.o.d -o CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.o -c /Users/nsdong2/Downloads/lemon-1.3.1/tools/dimacs-to-lgf.cc

tools/CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.i"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/tools && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nsdong2/Downloads/lemon-1.3.1/tools/dimacs-to-lgf.cc > CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.i

tools/CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.s"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/tools && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nsdong2/Downloads/lemon-1.3.1/tools/dimacs-to-lgf.cc -o CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.s

# Object files for target dimacs-to-lgf
dimacs__to__lgf_OBJECTS = \
"CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.o"

# External object files for target dimacs-to-lgf
dimacs__to__lgf_EXTERNAL_OBJECTS =

tools/dimacs-to-lgf: tools/CMakeFiles/dimacs-to-lgf.dir/dimacs-to-lgf.cc.o
tools/dimacs-to-lgf: tools/CMakeFiles/dimacs-to-lgf.dir/build.make
tools/dimacs-to-lgf: lemon/libemon.a
tools/dimacs-to-lgf: tools/CMakeFiles/dimacs-to-lgf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/nsdong2/Downloads/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable dimacs-to-lgf"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dimacs-to-lgf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tools/CMakeFiles/dimacs-to-lgf.dir/build: tools/dimacs-to-lgf
.PHONY : tools/CMakeFiles/dimacs-to-lgf.dir/build

tools/CMakeFiles/dimacs-to-lgf.dir/clean:
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/tools && $(CMAKE_COMMAND) -P CMakeFiles/dimacs-to-lgf.dir/cmake_clean.cmake
.PHONY : tools/CMakeFiles/dimacs-to-lgf.dir/clean

tools/CMakeFiles/dimacs-to-lgf.dir/depend:
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nsdong2/Downloads/lemon-1.3.1 /Users/nsdong2/Downloads/lemon-1.3.1/tools /Users/nsdong2/Downloads/lemon-1.3.1/build /Users/nsdong2/Downloads/lemon-1.3.1/build/tools /Users/nsdong2/Downloads/lemon-1.3.1/build/tools/CMakeFiles/dimacs-to-lgf.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : tools/CMakeFiles/dimacs-to-lgf.dir/depend

