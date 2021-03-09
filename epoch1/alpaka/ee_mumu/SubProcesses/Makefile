LIBDIR   = ../../lib
TOOLSDIR = ../../../../../tools
INCFLAGS = -I. -I../../src -I$(TOOLSDIR)
MODELLIB = model_sm
OPTFLAGS = -O3 -g
CXXFLAGS = $(OPTFLAGS) -std=c++14 $(INCFLAGS) -Wall -Wshadow -Wextra $(MGONGPU_CONFIG) -DFORCE_AOS
LIBFLAGS = -L$(LIBDIR) -l$(MODELLIB)
CXX     ?= g++

#Export CUDA_HOME to select a cuda installation
ifndef CUDA_HOME
  NVCC ?= $(shell which nvcc 2>/dev/null)
  ifneq ($(NVCC),)
    # NVCC is in the PATH or set explicitly
    CUDA_HOME  = $(patsubst %bin/nvcc,%,$(NVCC))
    CUDA_HOME := $(warning No CUDA_HOME exported. Using "$(CUDA_HOME)") $(CUDA_HOME)
  endif
endif

CUARCHNUM   = 70
#CUARCHNUM=61 #(For Pascal Architecture Cards)
USE_NVTX   ?= -DUSE_NVTX
CUARCHFLAGS = -arch=compute_$(CUARCHNUM)
CUINC       = -I$(CUDA_HOME)/include/
CULIBFLAGS  = -L$(CUDA_HOME)/lib64/ -lcuda
CUFLAGS     = $(OPTFLAGS) -std=c++14 $(INCFLAGS) $(CUINC) $(USE_NVTX) $(CUARCHFLAGS) -use_fast_math -lineinfo $(MGONGPU_CONFIG)


ifeq ($(shell test -f /opt/rocm/bin/hipcc; echo $$?),0)
  HIP_CC ?= /opt/rocm/bin/hipcc
endif
HIPFLAGS    = $(OPTFLAGS) -std=c++14 $(INCFLAGS) $(MGONGPU_CONFIG) -I/opt/rocm/hiprand/include -I/opt/rocm/rocrand/include

#Setup compiler
GCCV := $(shell g++ --version | head -1 | cut -d' ' -f3 | sed -e 's/\..*//')
ifeq ($(GCCV),4)
  $(error g++ compiler version 4 detected, which is too low)
endif

#Export ALPAKA_HOME to select an alpaka installation
ifndef ALPAKA_HOME
  ALPAKA_HOME := $(shell pwd)/alpaka
#  BOOSTINC ?= $(shell sed alpaka-repo/build/CMakeCache.txt -e '/^Boost_INCLUDE_DIR:PATH=/ !d; s/Boost_INCLUDE_DIR:PATH=//')
#  BOOSTINC := /usr/include/boost169
  BOOSTBASE := $(shell pwd)/boost
  BOOSTINC := $(shell pwd)/boost/include
endif

CMAKE3_BIN ?= $(shell which cmake3 2>/dev/null)
ifneq ($(CMAKE3_BIN),)
  CMAKE3VER := $(shell cmake3 --version | head -1 | sed -e 's/.*version //' -e 's/^3*\.\([0-9]*\).*/\1/')
  ifeq ($(shell test $(CMAKE3VER) -lt 15; echo $$?),0)
    CMAKE3_BIN :=
  endif
endif
ifeq ($(CMAKE3_BIN),)
  CMAKE_INSTALL := $(shell pwd)/cmake
  CMAKE3_BIN := $(CMAKE_INSTALL)/bin/cmake
  cmake3dep := cmake/bin/cmake
endif

#Export CUPLA_HOME to select the cupla source code (it will need to be build according to which backend we want)
ifndef CUPLA_HOME
  CUPLA_HOME := $(shell pwd)/cupla-repo
endif

ifeq ($(BOOSTINC),)
  $(error BOOSTINC not set)
endif

# build will need to contain various builds of cupla for each backend in which we're interested
CUPLA_BUILD_BASE := $(shell pwd)/cupla

# alpaka/cupla defined
ALPINC := $(ALPAKA_HOME)/usr/local/include
CUPLAINC := $(CUPLA_HOME)/include
CUPLASRC := $(CUPLA_HOME)/src
CUPLACPUFLAGS := -DALPAKA_ACC_CPU_B_OMP2_T_SEQ_ENABLED -fopenmp
CUPLACUDAFLAGS := -DALPAKA_ACC_GPU_CUDA_ENABLED
CUPLAHIPFLAGS := -DALPAKA_ACC_GPU_HIP_ENABLED

cxx_main=check_cpu.exe
cxx_cpu_objects=$(CUPLA_BUILD_BASE)/cpu/common.o $(CUPLA_BUILD_BASE)/cpu/device.o $(CUPLA_BUILD_BASE)/cpu/event.o $(CUPLA_BUILD_BASE)/cpu/memory.o $(CUPLA_BUILD_BASE)/cpu/stream.o $(CUPLA_BUILD_BASE)/cpu/manager/Driver.o

ifneq ($(NVCC),)
  cxx_main += check_gpu.exe
  cxx_gpu_objects = $(CUPLA_BUILD_BASE)/gpu/common.o $(CUPLA_BUILD_BASE)/gpu/device.o $(CUPLA_BUILD_BASE)/gpu/event.o $(CUPLA_BUILD_BASE)/gpu/memory.o $(CUPLA_BUILD_BASE)/gpu/stream.o $(CUPLA_BUILD_BASE)/gpu/manager/Driver.o
endif

ifneq ($(HIP_CC),)
  cxx_main += check_gpuhip.exe
  cxx_gpu_objects = $(CUPLA_BUILD_BASE)/gpuhip/common.o $(CUPLA_BUILD_BASE)/gpuhip/device.o $(CUPLA_BUILD_BASE)/gpuhip/event.o $(CUPLA_BUILD_BASE)/gpuhip/memory.o $(CUPLA_BUILD_BASE)/gpuhip/stream.o $(CUPLA_BUILD_BASE)/gpuhip/manager/Driver.o
endif

.PHONY: all debug
all: $(cmake3dep)
	if [ ! -d $(ALPAKA_HOME) ]; then \
		if [ ! -d $(BOOSTINC) ]; then \
			wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz; \
			tar zxf boost_1_75_0.tar.gz; \
			cd boost_1_75_0; \
			./bootstrap.sh; \
			./b2 install --prefix=$(BOOSTBASE); \
			cd ..; \
		fi; \
		if [ ! -d ./alpaka-repo ]; then \
			git clone -b 0.5.0 https://github.com/alpaka-group/alpaka.git alpaka-repo; \
			if [ $$? -ne 0 ]; then \
				echo "could not checkout alpaka"; \
				exit 1; \
			fi; \
		fi; \
		cd alpaka-repo; mkdir -p build; cd build; CMAKE_INCLUDE_PATH=$(BOOSTINC) $(CMAKE3_BIN) ..; make; make DESTDIR=$(ALPAKA_HOME) install; \
		sed -i 's/DimInt<3u>, uint32_t/DimInt<3u>, size_t/g' $(ALPAKA_HOME)/usr/local/include/alpaka/mem/buf/uniformCudaHip/Copy.hpp; \
	fi
	if [ ! -d $(CUPLA_HOME) ]; then \
		git clone -b 0.2.0 https://github.com/alpaka-group/cupla.git $(CUPLA_HOME); \
		if [ $$? -ne 0 ]; then \
			echo "could not checkout cupla"; \
			exit 1; \
		fi; \
	fi
#	cd ../../src && make
	mkdir -p cupla/cpu/manager cupla/gpu/manager cupla/gpuhip/manager
	make $(cxx_main)

debug: CXXFLAGS:=$(filter-out -O3,$(CXXFLAGS))
debug: CXXFLAGS += -O0 -DDEBUG2
debug: MAKEDEBUG := debug
debug: $(cxx_main)

$(LIBDIR)/lib$(MODELLIB).a:
	@cd ../../src && make ALPINC=$(ALPINC) BOOSTINC=$(BOOSTINC) $(MAKEDEBUG)

$(CUPLA_BUILD_BASE)/cpu/%.o: $(CUPLASRC)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(ALPINC) -I$(CUPLAINC) -I$(BOOSTINC) $(CUPLACPUFLAGS) -c $< -o $@

$(CUPLA_BUILD_BASE)/gpu/%.o: $(CUPLASRC)/%.cpp
	$(NVCC) -x cu $(CPPFLAGS) $(CUFLAGS) -I$(ALPINC) -I$(CUPLAINC) -I$(BOOSTINC) $(CUPLACUDAFLAGS) -c $< -o $@

$(CUPLA_BUILD_BASE)/gpuhip/%.o: $(CUPLASRC)/%.cpp
	$(HIP_CC) $(CPPFLAGS) $(HIPFLAGS) -I$(ALPINC) -I$(CUPLAINC) -I$(BOOSTINC) $(CUPLAHIPFLAGS) -c $< -o $@

check_cpu.o: check.cc CPPProcess.cc *.h ../../src/*.h ../../src/rambo.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(ALPINC) -I$(CUPLAINC) -I$(BOOSTINC) $(CUPLACPUFLAGS) -c check.cc -o $@

check_gpu.o: check.cc CPPProcess.cc *.h ../../src/*.h ../../src/rambo.cc
	$(NVCC) -x cu $(CPPFLAGS) $(CUFLAGS) -I$(ALPINC) -I$(CUPLAINC) -I$(BOOSTINC) $(CUPLACUDAFLAGS) -c check.cc -o $@

check_gpuhip.o: check.cc CPPProcess.cc *.h ../../src/*.h ../../src/rambo.cc
	$(HIP_CC) $(CPPFLAGS) $(HIPFLAGS) -I$(ALPINC) -I$(CUPLAINC) -I$(BOOSTINC) $(CUPLAHIPFLAGS) -c check.cc -o $@

check_cpu.exe: $(LIBDIR)/lib$(MODELLIB).a $(cxx_cpu_objects) check_cpu.o
	$(CXX) -o $@ $(cxx_cpu_objects) check_cpu.o -ldl $(LIBFLAGS) $(CUPLACPUFLAGS)

check_gpu.exe: $(LIBDIR)/lib$(MODELLIB).a $(cxx_gpu_objects) check_gpu.o
	$(NVCC) -o $@ $(cxx_gpu_objects) check_gpu.o -ldl $(LIBFLAGS) $(CULIBFLAGS) $(CUPLACUDAFLAGS)

check_gpuhip.exe: $(LIBDIR)/lib$(MODELLIB).a $(cxx_gpu_objects) check_gpuhip.o
	$(HIP_CC) -o $@ $(cxx_gpu_objects) check_gpuhip.o -ldl $(LIBFLAGS) $(CUPLAHIPFLAGS)

cmake/bin/cmake:
	wget https://github.com/Kitware/CMake/releases/download/v3.19.2/cmake-3.19.2.tar.gz
	tar zxf cmake-3.19.2.tar.gz
	cd cmake-3.19.2 && ./bootstrap --prefix=$(CMAKE_INSTALL) -- -DCMAKE_USE_OPENSSL:BOOL=OFF -DCMAKE_BUILD_TYPE:STRING=Release && make && make install

.PHONY: distclean clean perf test force
clean:
	cd ../../src && make clean
	rm -f *.o *.exe cupla/*/*.o cupla/*/manager/*.o

distclean:
	cd ../../src && make clean
	rm -f *.o *.exe cupla/*/*.o cupla/*/manager/*.o
	rm -rf alpaka-repo alpaka cupla-repo cupla cmake-3.19.2.tar.gz cmake-3.19.2 cmake boost boost_1_75_0 boost_1_75_0.tar.gz

perf: force
	make clean && make
	time ./check_gpu.exe -p 16384 32 12 && date

test: force
	./check_cpu.exe -v 32 1 1

force:
