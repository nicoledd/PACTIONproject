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
CMAKE_SOURCE_DIR = /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build

# Utility rule file for update-external-tags.

# Include any custom commands dependencies for this target.
include doc/CMakeFiles/update-external-tags.dir/compiler_depend.make

# Include the progress variables for this target.
include doc/CMakeFiles/update-external-tags.dir/progress.make

doc/CMakeFiles/update-external-tags:
	cd /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build/doc && /opt/homebrew/bin/wget -N http://gcc.gnu.org/onlinedocs/gcc-4.7.3/libstdc++/api/libstdc++.tag

update-external-tags: doc/CMakeFiles/update-external-tags
update-external-tags: doc/CMakeFiles/update-external-tags.dir/build.make
.PHONY : update-external-tags

# Rule to build all files generated by this target.
doc/CMakeFiles/update-external-tags.dir/build: update-external-tags
.PHONY : doc/CMakeFiles/update-external-tags.dir/build

doc/CMakeFiles/update-external-tags.dir/clean:
	cd /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build/doc && $(CMAKE_COMMAND) -P CMakeFiles/update-external-tags.dir/cmake_clean.cmake
.PHONY : doc/CMakeFiles/update-external-tags.dir/clean

doc/CMakeFiles/update-external-tags.dir/depend:
	cd /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1 /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/doc /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build/doc /Users/nsdong2/Documents/PACTIONproject/lemon-1.3.1/build/doc/CMakeFiles/update-external-tags.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : doc/CMakeFiles/update-external-tags.dir/depend

