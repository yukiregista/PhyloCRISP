# Version of booster
GIT_VERSION := $(shell git describe --abbrev=10 --dirty --always --tags)

UNAME := $(shell uname)
CC_VERSION = $(shell $(CC) --version 2>/dev/null | head -n 1)

#Use -DDEBUG to print debug output.  Use -DNDEBUG to ignore asserts.

# FOR DEBUGGING
# CFLAGS = -Wall -g -DVERSION=\"$(GIT_VERSION)\" -DDEBUG -fno-stack-protector
# CFLAGS_OMP = -Wall -g -fopenmp -DDEBUG


OMP_CFLAGS = -fopenmp
OMP_LIBS =
CFLAGS = -Wall -g -DVERSION=\"$(GIT_VERSION)\" -DDEBUG
CFLAGS_OMP = -Wall -g $(OMP_CFLAGS) -DDEBUG

# assertはされるがDB_TRACEなどのプリントはされないバージョン
CFLAGS = -Wall -g -DVERSION=\"$(GIT_VERSION)\"
CFLAGS_OMP = -Wall -g $(OMP_CFLAGS)

# 本番
# CFLAGS = -Wall -O3 -DVERSION=\"$(GIT_VERSION)\" -DNDEBUG
# CFLAGS_OMP = -Wall -O3 -fopenmp -DNDEBUG
# CFLAGS_SHARED = -Wall -O3 -fopenmp -fPIC -DVERSION=\"$(GIT_VERSION)\" -DNDEBUG  # Added
# CFLAGS_SHARED = -Wall -O3 -fopenmp -fPIC -DVERSION=\"$(GIT_VERSION)\" -DNDEBUG  # Added




# Compiler: gcc
ifeq ($(cross),win32)
        CC = i686-w64-mingw32-gcc
else
	ifeq ($(cross),win64)
		CC = x86_64-w64-mingw32-gcc
	else
		ifeq ($(cross),linux32)
			CFLAGS_OMP += -m32
			CFLAGS += -m32
			CFLAGS_SHARED += -m32  # Added
		else
			ifeq ($(UNAME),Darwin)
				CC = clang
			else
				CC = gcc
			endif
		endif
	endif
endif

ifeq ($(UNAME),Darwin)
	# Apple Clang does not support plain '-fopenmp' or '-static-libgcc'.
	# Use Homebrew libomp flags when compiling with clang.
	ifneq (,$(findstring clang,$(CC_VERSION)))
		LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
		ifneq ($(LIBOMP_PREFIX),)
			OMP_CFLAGS = -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
			OMP_LIBS = -L$(LIBOMP_PREFIX)/lib -lomp
			CFLAGS_OMP = -Wall -g $(OMP_CFLAGS)
		else
			$(error libomp not found. Run 'brew install libomp' or build with Homebrew GCC, e.g. 'make CC=gcc-14')
		endif
	endif
endif

LIBS = -lm
OBJS = hashtables_bfields.o  tree.o stats.o prng.o hashmap.o version.o sort.o io.o tree_utils.o bitset_index.o rapid_transfer.o heavy_paths.o debug.o kludge.o node_stack.o split.o M1M2_hashmap.o transfer_distance.o

# default target
ALL = booster

INSTALL_PATH=$$HOME/bin/

all : $(ALL)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

# ****
# the "booster" supports. Needs ref tree and bt trees.
# ****
booster: $(OBJS) booster.c
	$(CC) $(CFLAGS_OMP) -o $@ $^ $(LIBS) $(OMP_LIBS)

# Build the shared library (new target)
shared: libbooster.so  # Added

libbooster.so: $(OBJS) booster.o  # Added
	$(CC) $(CFLAGS_SHARED) -shared -o $@ $^ $(LIBS)



# ****
# TESTS
# ****
tests: $(OBJS) test.c
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

test : tests
	./tests

.PHONY: clean

clean:
	rm -f *~ *.o $(ALL) tests 
	rm -rf *.dSYM

install: all
	mkdir -p $(INSTALL_PATH)
	cp $(ALL) $(INSTALL_PATH)

uninstall:
	rm $(addprefix $(INSTALL_PATH),$(ALL))
