# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.25.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.25.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/phamtrung/Desktop/Code/Mpboot1/mpboot

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/phamtrung/Desktop/Code/Mpboot1/mpboot

# Include any dependencies generated for this target.
include sprng/CMakeFiles/sprng.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include sprng/CMakeFiles/sprng.dir/compiler_depend.make

# Include the progress variables for this target.
include sprng/CMakeFiles/sprng.dir/progress.make

# Include the compile flags for this target's objects.
include sprng/CMakeFiles/sprng.dir/flags.make

sprng/CMakeFiles/sprng.dir/lcg64.c.o: sprng/CMakeFiles/sprng.dir/flags.make
sprng/CMakeFiles/sprng.dir/lcg64.c.o: sprng/lcg64.c
sprng/CMakeFiles/sprng.dir/lcg64.c.o: sprng/CMakeFiles/sprng.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot1/mpboot/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object sprng/CMakeFiles/sprng.dir/lcg64.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT sprng/CMakeFiles/sprng.dir/lcg64.c.o -MF CMakeFiles/sprng.dir/lcg64.c.o.d -o CMakeFiles/sprng.dir/lcg64.c.o -c /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/lcg64.c

sprng/CMakeFiles/sprng.dir/lcg64.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sprng.dir/lcg64.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/lcg64.c > CMakeFiles/sprng.dir/lcg64.c.i

sprng/CMakeFiles/sprng.dir/lcg64.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sprng.dir/lcg64.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/lcg64.c -o CMakeFiles/sprng.dir/lcg64.c.s

sprng/CMakeFiles/sprng.dir/makeseed.c.o: sprng/CMakeFiles/sprng.dir/flags.make
sprng/CMakeFiles/sprng.dir/makeseed.c.o: sprng/makeseed.c
sprng/CMakeFiles/sprng.dir/makeseed.c.o: sprng/CMakeFiles/sprng.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot1/mpboot/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object sprng/CMakeFiles/sprng.dir/makeseed.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT sprng/CMakeFiles/sprng.dir/makeseed.c.o -MF CMakeFiles/sprng.dir/makeseed.c.o.d -o CMakeFiles/sprng.dir/makeseed.c.o -c /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/makeseed.c

sprng/CMakeFiles/sprng.dir/makeseed.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sprng.dir/makeseed.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/makeseed.c > CMakeFiles/sprng.dir/makeseed.c.i

sprng/CMakeFiles/sprng.dir/makeseed.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sprng.dir/makeseed.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/makeseed.c -o CMakeFiles/sprng.dir/makeseed.c.s

sprng/CMakeFiles/sprng.dir/memory.c.o: sprng/CMakeFiles/sprng.dir/flags.make
sprng/CMakeFiles/sprng.dir/memory.c.o: sprng/memory.c
sprng/CMakeFiles/sprng.dir/memory.c.o: sprng/CMakeFiles/sprng.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot1/mpboot/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object sprng/CMakeFiles/sprng.dir/memory.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT sprng/CMakeFiles/sprng.dir/memory.c.o -MF CMakeFiles/sprng.dir/memory.c.o.d -o CMakeFiles/sprng.dir/memory.c.o -c /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/memory.c

sprng/CMakeFiles/sprng.dir/memory.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sprng.dir/memory.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/memory.c > CMakeFiles/sprng.dir/memory.c.i

sprng/CMakeFiles/sprng.dir/memory.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sprng.dir/memory.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/memory.c -o CMakeFiles/sprng.dir/memory.c.s

sprng/CMakeFiles/sprng.dir/store.c.o: sprng/CMakeFiles/sprng.dir/flags.make
sprng/CMakeFiles/sprng.dir/store.c.o: sprng/store.c
sprng/CMakeFiles/sprng.dir/store.c.o: sprng/CMakeFiles/sprng.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot1/mpboot/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object sprng/CMakeFiles/sprng.dir/store.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT sprng/CMakeFiles/sprng.dir/store.c.o -MF CMakeFiles/sprng.dir/store.c.o.d -o CMakeFiles/sprng.dir/store.c.o -c /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/store.c

sprng/CMakeFiles/sprng.dir/store.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sprng.dir/store.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/store.c > CMakeFiles/sprng.dir/store.c.i

sprng/CMakeFiles/sprng.dir/store.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sprng.dir/store.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/store.c -o CMakeFiles/sprng.dir/store.c.s

sprng/CMakeFiles/sprng.dir/primes-lcg64.c.o: sprng/CMakeFiles/sprng.dir/flags.make
sprng/CMakeFiles/sprng.dir/primes-lcg64.c.o: sprng/primes-lcg64.c
sprng/CMakeFiles/sprng.dir/primes-lcg64.c.o: sprng/CMakeFiles/sprng.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot1/mpboot/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object sprng/CMakeFiles/sprng.dir/primes-lcg64.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT sprng/CMakeFiles/sprng.dir/primes-lcg64.c.o -MF CMakeFiles/sprng.dir/primes-lcg64.c.o.d -o CMakeFiles/sprng.dir/primes-lcg64.c.o -c /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/primes-lcg64.c

sprng/CMakeFiles/sprng.dir/primes-lcg64.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sprng.dir/primes-lcg64.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/primes-lcg64.c > CMakeFiles/sprng.dir/primes-lcg64.c.i

sprng/CMakeFiles/sprng.dir/primes-lcg64.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sprng.dir/primes-lcg64.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/primes-lcg64.c -o CMakeFiles/sprng.dir/primes-lcg64.c.s

sprng/CMakeFiles/sprng.dir/checkid.c.o: sprng/CMakeFiles/sprng.dir/flags.make
sprng/CMakeFiles/sprng.dir/checkid.c.o: sprng/checkid.c
sprng/CMakeFiles/sprng.dir/checkid.c.o: sprng/CMakeFiles/sprng.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot1/mpboot/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object sprng/CMakeFiles/sprng.dir/checkid.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT sprng/CMakeFiles/sprng.dir/checkid.c.o -MF CMakeFiles/sprng.dir/checkid.c.o.d -o CMakeFiles/sprng.dir/checkid.c.o -c /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/checkid.c

sprng/CMakeFiles/sprng.dir/checkid.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/sprng.dir/checkid.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/checkid.c > CMakeFiles/sprng.dir/checkid.c.i

sprng/CMakeFiles/sprng.dir/checkid.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/sprng.dir/checkid.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/checkid.c -o CMakeFiles/sprng.dir/checkid.c.s

# Object files for target sprng
sprng_OBJECTS = \
"CMakeFiles/sprng.dir/lcg64.c.o" \
"CMakeFiles/sprng.dir/makeseed.c.o" \
"CMakeFiles/sprng.dir/memory.c.o" \
"CMakeFiles/sprng.dir/store.c.o" \
"CMakeFiles/sprng.dir/primes-lcg64.c.o" \
"CMakeFiles/sprng.dir/checkid.c.o"

# External object files for target sprng
sprng_EXTERNAL_OBJECTS =

sprng/libsprng.a: sprng/CMakeFiles/sprng.dir/lcg64.c.o
sprng/libsprng.a: sprng/CMakeFiles/sprng.dir/makeseed.c.o
sprng/libsprng.a: sprng/CMakeFiles/sprng.dir/memory.c.o
sprng/libsprng.a: sprng/CMakeFiles/sprng.dir/store.c.o
sprng/libsprng.a: sprng/CMakeFiles/sprng.dir/primes-lcg64.c.o
sprng/libsprng.a: sprng/CMakeFiles/sprng.dir/checkid.c.o
sprng/libsprng.a: sprng/CMakeFiles/sprng.dir/build.make
sprng/libsprng.a: sprng/CMakeFiles/sprng.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot1/mpboot/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking C static library libsprng.a"
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && $(CMAKE_COMMAND) -P CMakeFiles/sprng.dir/cmake_clean_target.cmake
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sprng.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
sprng/CMakeFiles/sprng.dir/build: sprng/libsprng.a
.PHONY : sprng/CMakeFiles/sprng.dir/build

sprng/CMakeFiles/sprng.dir/clean:
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng && $(CMAKE_COMMAND) -P CMakeFiles/sprng.dir/cmake_clean.cmake
.PHONY : sprng/CMakeFiles/sprng.dir/clean

sprng/CMakeFiles/sprng.dir/depend:
	cd /Users/phamtrung/Desktop/Code/Mpboot1/mpboot && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/phamtrung/Desktop/Code/Mpboot1/mpboot /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng /Users/phamtrung/Desktop/Code/Mpboot1/mpboot /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng /Users/phamtrung/Desktop/Code/Mpboot1/mpboot/sprng/CMakeFiles/sprng.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : sprng/CMakeFiles/sprng.dir/depend
