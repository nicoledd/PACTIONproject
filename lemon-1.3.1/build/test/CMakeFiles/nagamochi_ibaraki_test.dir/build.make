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
include test/CMakeFiles/nagamochi_ibaraki_test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/nagamochi_ibaraki_test.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/nagamochi_ibaraki_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/nagamochi_ibaraki_test.dir/flags.make

test/CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.o: test/CMakeFiles/nagamochi_ibaraki_test.dir/flags.make
test/CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.o: /Users/nsdong2/Downloads/lemon-1.3.1/test/nagamochi_ibaraki_test.cc
test/CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.o: test/CMakeFiles/nagamochi_ibaraki_test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/nsdong2/Downloads/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.o"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.o -MF CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.o.d -o CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.o -c /Users/nsdong2/Downloads/lemon-1.3.1/test/nagamochi_ibaraki_test.cc

test/CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.i"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nsdong2/Downloads/lemon-1.3.1/test/nagamochi_ibaraki_test.cc > CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.i

test/CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.s"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nsdong2/Downloads/lemon-1.3.1/test/nagamochi_ibaraki_test.cc -o CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.s

# Object files for target nagamochi_ibaraki_test
nagamochi_ibaraki_test_OBJECTS = \
"CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.o"

# External object files for target nagamochi_ibaraki_test
nagamochi_ibaraki_test_EXTERNAL_OBJECTS =

test/nagamochi_ibaraki_test: test/CMakeFiles/nagamochi_ibaraki_test.dir/nagamochi_ibaraki_test.cc.o
test/nagamochi_ibaraki_test: test/CMakeFiles/nagamochi_ibaraki_test.dir/build.make
test/nagamochi_ibaraki_test: lemon/libemon.a
test/nagamochi_ibaraki_test: test/CMakeFiles/nagamochi_ibaraki_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/nsdong2/Downloads/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable nagamochi_ibaraki_test"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nagamochi_ibaraki_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/nagamochi_ibaraki_test.dir/build: test/nagamochi_ibaraki_test
.PHONY : test/CMakeFiles/nagamochi_ibaraki_test.dir/build

test/CMakeFiles/nagamochi_ibaraki_test.dir/clean:
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -P CMakeFiles/nagamochi_ibaraki_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/nagamochi_ibaraki_test.dir/clean

test/CMakeFiles/nagamochi_ibaraki_test.dir/depend:
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nsdong2/Downloads/lemon-1.3.1 /Users/nsdong2/Downloads/lemon-1.3.1/test /Users/nsdong2/Downloads/lemon-1.3.1/build /Users/nsdong2/Downloads/lemon-1.3.1/build/test /Users/nsdong2/Downloads/lemon-1.3.1/build/test/CMakeFiles/nagamochi_ibaraki_test.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/CMakeFiles/nagamochi_ibaraki_test.dir/depend

