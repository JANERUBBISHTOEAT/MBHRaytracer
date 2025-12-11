JLEVEL=-j4
all: debug release openmp

# Auto-detect BOOST_ROOT: prefer sysroot if exists, otherwise use system boost
BOOST_ROOT ?= $(shell if [ -d "${HOME}/sysroot" ]; then echo "${HOME}/sysroot"; else echo "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc11/boost/1.80.0"; fi)
CARGS=-DBOOST_ROOT=$(BOOST_ROOT)  -DCMAKE_MODULE_PATH=$(BOOST_ROOT)/lib64/cmake/fmt/

### debug: Generate debug version of the code
debug: dirs
	cd debug && cmake ${CARGS} .. && make $(JLEVEL)

### openmp: Generate OMP shared-memory parallelism version of the code
openmp: dirs
	cd openmp && cmake ${CARGS} -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=False -DENABLE_OPENMP=True .. && make $(JLEVEL)

### release: Generate MPI+OMP parallelism version of the code
release: dirs
	cd release && cmake ${CARGS} -DCMAKE_BUILD_TYPE=Release -DENABLE_MPI=True -DENABLE_OPENMP=True .. && make $(JLEVEL)

dirs:
	mkdir -p debug release openmp


### docs: Generate documentation with Doxygen
docs:
	doxygen Doxyfile

### help: Makefile's descriptive rules help
help: Makefile
	@sed -n 's/^###//p' $<

### clean: Remove auxiliary object files
clean:
	cd debug && make clean
	cd release && make clean
	cd openmp && make clean

### distclean: Remove all builds
distclean:
	-rm -rf debug release openmp

#cppcheck:
#	cppcheck --language=c++ --std=c++17 main.cc mov_avg.cc mov_avg.h
format:
	find src -iname '*.h' -o -iname '*.cc' | xargs clang-format -i
