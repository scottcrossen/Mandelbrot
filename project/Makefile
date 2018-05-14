# Makefile by Scott Leland Crossen
# @Copyright 2017 Scott Leland Crossen

# ==== Definitions ====

BYUNETID = slxn42
TMPDIR = tmp
SOURCEDIR = src
RELDIR = rel
COMP = gcc
BINOBJS = ${TMPDIR}/mandelbrot.o
MAIN = ${RELDIR}/mandelbrot
COMPFLAGS = -g -fopenmp -Wall
LINKFLAGS =

# ==== Recipes ====

.PHONY: default
.DEFAULT: default
default: all

.PHONY: all
all: $(MAIN)

.PHONY: clean
clean:
	rm -rf $(TMPDIR) $(RELDIR) $(MAIN)

# ==== Targets ====

.PRECIOUS: ${TMPDIR}/mandelbrot.o
${TMPDIR}/mandelbrot.o: ${SOURCEDIR}/mandelbrot.c
	mkdir -p $(@D)
	$(COMP) -c $< -o $@ $(COMPFLAGS)

.PRECIOUS: ${TMPDIR}/%.o
${TMPDIR}/%.o: ${SOURCEDIR}/%.c ${SOURCEDIR}/%.h
	mkdir -p $(@D)
	$(COMP) -c $< -o $@ $(COMPFLAGS)

.PRECIOUS: ${RELDIR}/%
${RELDIR}/%: $(BINOBJS)
	mkdir -p $(@D)
	$(COMP) $^ -o $@ $(COMPFLAGS) $(LINKFLAGS)
