FC = mpif90
FFLAGS = -O2
# Homebrew gfortran on macOS needs the SDK to link libSystem
SDKROOT ?= $(shell xcrun --show-sdk-path 2>/dev/null)

.PHONY: all clean

all: wham1d

wham1d: wham_1D.f90
	SDKROOT=$(SDKROOT) $(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f wham1d *.o *.mod
