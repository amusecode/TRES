ifneq (,$(filter-out clean distclean, $(MAKECMDGOALS)))
include support/config.mk
endif


# Detecting dependencies
support/config.mk:
	cd support && ./configure


# Building the code into a static library
##### Remove anything not needed #####

CFLAGS += 
LDFLAGS += 

OBJS = src/main_code.o src/ODE_system.o src/helper_routines.o src/tidal_friction_parameters.o src/cvode/cvode.o src/cvode/cvode_dense.o src/cvode/cvode_direct.o src/cvode/cvode_io.o src/cvode/cvode_diag.o src/cvode/cvode_spils.o src/cvode/cvode_spgmr.o src/cvode/cvode_spbcgs.o src/cvode/cvode_sptfqmr.o src/cvode/nvector_serial.o src/cvode/sundials_dense.o src/cvode/sundials_direct.o src/cvode/sundials_math.o src/cvode/sundials_nvector.o src/cvode/sundials_spgmr.o src/cvode/sundials_iterative.o src/cvode/sundials_spbcgs.o src/cvode/sundials_sptfqmr.o

# Building the workers
seculartriple_worker.h: interface.py
	amusifier --type=h interface.py SeculartripleInterface -o $@

seculartriple_worker.cc: interface.py
	amusifier --type=c interface.py SeculartripleInterface -o $@

seculartriple_worker.o: seculartriple_worker.cc
	$(MPICXX) -Isrc -c -o $@ $(CFLAGS) $<

seculartriple_worker: seculartriple_worker.o $(OBJS)
	$(MPICXX) -o $@ $(LDFLAGS) $^ $(LDLIBS)




# Which packages contain which workers?
amuse-seculartriple_contains: seculartriple_worker


# Building and installing packages
develop-%: %_contains
	support/shared/uninstall.sh $*
	python -m pip install -e packages/$*

install-%: %_contains
	support/shared/uninstall.sh $*
	python -m pip install packages/$*

package-%: %_contains
	python3 -m pip install -vv --no-cache-dir --no-deps --no-build-isolation --prefix ${PREFIX} packages/$*

test-%:
	cd packages/$* && pytest


# Cleaning up
.PHONY: clean
clean:
	$(MAKE) -C src clean
	rm -rf *.o *.mod *worker*

.PHONY: distclean
distclean: clean
	rm -f support/config.mk support/config.log support/config.status
	rm -rf support/autom4te.cache

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.c.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $<
