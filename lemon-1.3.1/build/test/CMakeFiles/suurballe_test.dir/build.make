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
include test/CMakeFiles/suurballe_test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/suurballe_test.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/suurballe_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/suurballe_test.dir/flags.make

test/CMakeFiles/suurballe_test.dir/suurballe_test.cc.o: test/CMakeFiles/suurballe_test.dir/flags.make
test/CMakeFiles/suurballe_test.dir/suurballe_test.cc.o: /Users/nsdong2/Downloads/lemon-1.3.1/test/suurballe_test.cc
test/CMakeFiles/suurballe_test.dir/suurballe_test.cc.o: test/CMakeFiles/suurballe_test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/nsdong2/Downloads/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/suurballe_test.dir/suurballe_test.cc.o"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/suurballe_test.dir/suurballe_test.cc.o -MF CMakeFiles/suurballe_test.dir/suurballe_test.cc.o.d -o CMakeFiles/suurballe_test.dir/suurballe_test.cc.o -c /Users/nsdong2/Downloads/lemon-1.3.1/test/suurballe_test.cc

test/CMakeFiles/suurballe_test.dir/suurballe_test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/suurballe_test.dir/suurballe_test.cc.i"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/nsdong2/Downloads/lemon-1.3.1/test/suurballe_test.cc > CMakeFiles/suurballe_test.dir/suurballe_test.cc.i

test/CMakeFiles/suurballe_test.dir/suurballe_test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/suurballe_test.dir/suurballe_test.cc.s"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/nsdong2/Downloads/lemon-1.3.1/test/suurballe_test.cc -o CMakeFiles/suurballe_test.dir/suurballe_test.cc.s

# Object files for target suurballe_test
suurballe_test_OBJECTS = \
"CMakeFiles/suurballe_test.dir/suurballe_test.cc.o"

# External object files for target suurballe_test
suurballe_test_EXTERNAL_OBJECTS =

test/suurballe_test: test/CMakeFiles/suurballe_test.dir/suurballe_test.cc.o
test/suurballe_test: test/CMakeFiles/suurballe_test.dir/build.make
test/suurballe_test: lemon/libemon.a
test/suurballe_test: test/CMakeFiles/suurballe_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/nsdong2/Downloads/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable suurballe_test"
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/suurballe_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/suurballe_test.dir/build: test/suurballe_test
.PHONY : test/CMakeFiles/suurballe_test.dir/build

test/CMakeFiles/suurballe_test.dir/clean:
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -P CMakeFiles/suurballe_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/suurballe_test.dir/clean

test/CMakeFiles/suurballe_test.dir/depend:
	cd /Users/nsdong2/Downloads/lemon-1.3.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nsdong2/Downloads/lemon-1.3.1 /Users/nsdong2/Downloads/lemon-1.3.1/test /Users/nsdong2/Downloads/lemon-1.3.1/build /Users/nsdong2/Downloads/lemon-1.3.1/build/test /Users/nsdong2/Downloads/lemon-1.3.1/build/test/CMakeFiles/suurballe_test.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : test/CMakeFiles/suurballe_test.dir/depend

