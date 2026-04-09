# Makefile for Beam Section Tool
# Authors: Bruno Zilli & DeepSeek
# Licence: MIT

FC = gfortran
FFLAGS = -Wall -O2 -ffixed-form -ffixed-line-length-132
TARGETS = test_torsion

OBJS = src/compute_section_properties.o \
       src/read_section_mesh_unv.o \
       src/section_database.o \
       src/mesh_checker.o \
       src/torsion_j.o

all: $(TARGETS)

test_torsion: test/test_torsion.f $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ -llapack -lblas

src/%.o: src/%.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f src/*.o src/*.mod test_torsion

.PHONY: all clean
