# Makefile for ray3d and Shibutani code (theo)
include ../Makefile.compiler

CC  = /usr/bin/gcc

# Change this if the ray3d code is elsewhere:
RAY3D_DIR = .
THEO_DIR = .

# Change this if the executive is to place elsewhere
BIN_DIR = ./bin

default: all

all: ray3d theo four qlayer

ray3d:
	$(FCRS) -c $(RAY3D_DIR)/ray3d.f

theo:
	$(FCRS) -c $(THEO_DIR)/theo.f90
four:
	$(FCRS) -c $(THEO_DIR)/four1.f
qlayer:
	$(FCRS) -c $(THEO_DIR)/qlayer.f90

clean:
	rm -f *.o core
