# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5


CMakeFiles/coalesce.dir/constants.mod.proxy: CMakeFiles/coalesce.dir/src/constants.f90.o.provides
CMakeFiles/coalesce.dir/src/constants.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod mod/constants CMakeFiles/coalesce.dir/constants.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/coalesce.dir/src/constants.f90.o.provides.build
CMakeFiles/coalesce.dir/build: CMakeFiles/coalesce.dir/src/constants.f90.o.provides.build

CMakeFiles/coalesce.dir/src/globals.f90.o.requires: CMakeFiles/coalesce.dir/constants.mod.proxy
CMakeFiles/coalesce.dir/src/globals.f90.o: CMakeFiles/coalesce.dir/constants.mod.stamp
CMakeFiles/coalesce.dir/globals.mod.proxy: CMakeFiles/coalesce.dir/src/globals.f90.o.provides
CMakeFiles/coalesce.dir/src/globals.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod mod/globals CMakeFiles/coalesce.dir/globals.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/coalesce.dir/src/globals.f90.o.provides.build
CMakeFiles/coalesce.dir/build: CMakeFiles/coalesce.dir/src/globals.f90.o.provides.build

CMakeFiles/coalesce.dir/src/in_out.f90.o.requires: CMakeFiles/coalesce.dir/constants.mod.proxy
CMakeFiles/coalesce.dir/src/in_out.f90.o: CMakeFiles/coalesce.dir/constants.mod.stamp
CMakeFiles/coalesce.dir/src/in_out.f90.o: /home/kvetinka/lib/cmake/fson/../../../include/fson/fson.mod
CMakeFiles/coalesce.dir/src/in_out.f90.o: /home/kvetinka/lib/cmake/fson/../../../include/fson/fson_value_m.mod
CMakeFiles/coalesce.dir/src/in_out.f90.o.requires: CMakeFiles/coalesce.dir/globals.mod.proxy
CMakeFiles/coalesce.dir/src/in_out.f90.o: CMakeFiles/coalesce.dir/globals.mod.stamp
CMakeFiles/coalesce.dir/src/in_out.f90.o.requires: CMakeFiles/coalesce.dir/ioutils.mod.proxy
CMakeFiles/coalesce.dir/src/in_out.f90.o: CMakeFiles/coalesce.dir/ioutils.mod.stamp
CMakeFiles/coalesce.dir/in_out.mod.proxy: CMakeFiles/coalesce.dir/src/in_out.f90.o.provides
CMakeFiles/coalesce.dir/src/in_out.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod mod/in_out CMakeFiles/coalesce.dir/in_out.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/coalesce.dir/src/in_out.f90.o.provides.build
CMakeFiles/coalesce.dir/build: CMakeFiles/coalesce.dir/src/in_out.f90.o.provides.build

CMakeFiles/coalesce.dir/src/integration.f90.o.requires: CMakeFiles/coalesce.dir/constants.mod.proxy
CMakeFiles/coalesce.dir/src/integration.f90.o: CMakeFiles/coalesce.dir/constants.mod.stamp
CMakeFiles/coalesce.dir/src/integration.f90.o.requires: CMakeFiles/coalesce.dir/globals.mod.proxy
CMakeFiles/coalesce.dir/src/integration.f90.o: CMakeFiles/coalesce.dir/globals.mod.stamp
CMakeFiles/coalesce.dir/src/integration.f90.o.requires: CMakeFiles/coalesce.dir/in_out.mod.proxy
CMakeFiles/coalesce.dir/src/integration.f90.o: CMakeFiles/coalesce.dir/in_out.mod.stamp
CMakeFiles/coalesce.dir/src/integration.f90.o.requires: CMakeFiles/coalesce.dir/list.mod.proxy
CMakeFiles/coalesce.dir/src/integration.f90.o: CMakeFiles/coalesce.dir/list.mod.stamp
CMakeFiles/coalesce.dir/src/integration.f90.o.requires: CMakeFiles/coalesce.dir/model.mod.proxy
CMakeFiles/coalesce.dir/src/integration.f90.o: CMakeFiles/coalesce.dir/model.mod.stamp
CMakeFiles/coalesce.dir/integration.mod.proxy: CMakeFiles/coalesce.dir/src/integration.f90.o.provides
CMakeFiles/coalesce.dir/src/integration.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod mod/integration CMakeFiles/coalesce.dir/integration.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/coalesce.dir/src/integration.f90.o.provides.build
CMakeFiles/coalesce.dir/build: CMakeFiles/coalesce.dir/src/integration.f90.o.provides.build

CMakeFiles/coalesce.dir/ioutils.mod.proxy: CMakeFiles/coalesce.dir/src/ioutils.f90.o.provides
CMakeFiles/coalesce.dir/src/ioutils.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod mod/ioutils CMakeFiles/coalesce.dir/ioutils.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/coalesce.dir/src/ioutils.f90.o.provides.build
CMakeFiles/coalesce.dir/build: CMakeFiles/coalesce.dir/src/ioutils.f90.o.provides.build

CMakeFiles/coalesce.dir/list.mod.proxy: CMakeFiles/coalesce.dir/src/list.f90.o.provides
CMakeFiles/coalesce.dir/src/list.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod mod/list CMakeFiles/coalesce.dir/list.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/coalesce.dir/src/list.f90.o.provides.build
CMakeFiles/coalesce.dir/build: CMakeFiles/coalesce.dir/src/list.f90.o.provides.build

CMakeFiles/coalesce.dir/src/main.f90.o.requires: CMakeFiles/coalesce.dir/in_out.mod.proxy
CMakeFiles/coalesce.dir/src/main.f90.o: CMakeFiles/coalesce.dir/in_out.mod.stamp
CMakeFiles/coalesce.dir/src/main.f90.o.requires: CMakeFiles/coalesce.dir/integration.mod.proxy
CMakeFiles/coalesce.dir/src/main.f90.o: CMakeFiles/coalesce.dir/integration.mod.stamp

CMakeFiles/coalesce.dir/src/model.f90.o.requires: CMakeFiles/coalesce.dir/constants.mod.proxy
CMakeFiles/coalesce.dir/src/model.f90.o: CMakeFiles/coalesce.dir/constants.mod.stamp
CMakeFiles/coalesce.dir/src/model.f90.o.requires: CMakeFiles/coalesce.dir/globals.mod.proxy
CMakeFiles/coalesce.dir/src/model.f90.o: CMakeFiles/coalesce.dir/globals.mod.stamp
CMakeFiles/coalesce.dir/model.mod.proxy: CMakeFiles/coalesce.dir/src/model.f90.o.provides
CMakeFiles/coalesce.dir/src/model.f90.o.provides.build:
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod mod/model CMakeFiles/coalesce.dir/model.mod.stamp GNU
	$(CMAKE_COMMAND) -E touch CMakeFiles/coalesce.dir/src/model.f90.o.provides.build
CMakeFiles/coalesce.dir/build: CMakeFiles/coalesce.dir/src/model.f90.o.provides.build



