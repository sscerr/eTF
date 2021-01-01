include ../Machine

COMM= alltoall.o bcast.o cp_mype.o gather.o globalindex.o griddim.o \
 hangup.o localdim.o localindex.o ownerof.o reduce.o \
 scatter.o sendrecv.o startup.o syncronize.o

all: $(COMM)
	ar rv comm.a *.o


include ../Rules
