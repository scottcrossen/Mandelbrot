# Makefile by Scott Leland Crossen
# @Copyright 2017 Scott Leland Crossen

# ==== Definitions ====

BYUNETID = slxn42
BINDIR = bin
SOURCEDIR = src
RELDIR = rel
COMP = gcc
BINOBJS = ${BINDIR}/mandelbrot.o
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
	rm -rf $(BINDIR) $(RELDIR) $(MAIN)

# ==== Targets ====

.PRECIOUS: ${BINDIR}/mandelbrot.o
${BINDIR}/mandelbrot.o: ${SOURCEDIR}/mandelbrot.c
	mkdir -p $(@D)
	$(COMP) -c $< -o $@ $(COMPFLAGS)

.PRECIOUS: ${BINDIR}/%.o
${BINDIR}/%.o: ${SOURCEDIR}/%.c ${SOURCEDIR}/%.h
	mkdir -p $(@D)
	$(COMP) -c $< -o $@ $(COMPFLAGS)

.PRECIOUS: ${RELDIR}/%
${RELDIR}/%: $(BINOBJS)
	mkdir -p $(@D)
	$(COMP) $^ -o $@ $(COMPFLAGS) $(LINKFLAGS)
