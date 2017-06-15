CPP = g++
PROJ = CartPend
DEPS = $(PROJ).cpp
USERDIR = ./user

# Require SAC_PATH environment variable to lib dir
LIBDIR = /home/alex/Documents/code/C/SAC_v2/lib/src

# Prefer static libs to dynamic
LIBS = -Wl,-Bstatic \
	-lboost_system -lboost_timer -lboost_chrono -lrt -Wl,-Bdynamic
CFLAGS = -I/usr/include/eigen3/ -I$(USERDIR) -I$(LIBDIR)
CFLAGS += -std=c++11

# Default code optimization level
OPT = -O2

all: $(PROJ)

# Debug build with extra warning turned on and optimization off
debug: OPT = -O0 -g -m64
debug: CFLAGS += -Wall -Wextra -std=c++11
debug: $(PROJ)

$(PROJ): $(DEPS)
	$(CPP) $(CFLAGS) $^ -o $@ $(OPT) $(LIBS)

clean:
	rm -f *.o $(PROJ)
