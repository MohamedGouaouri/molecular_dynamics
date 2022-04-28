LIB        = -L. pthreads/initVeloFun.cpp
INCLUDE    = -I. -pthread
CFLAGS     = -O3 
EXEC       = MD
BUILD_PATH = build
CXX        = g++

${EXEC}: MD.cpp
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} MD.cpp -o ${BUILD_PATH}/${EXEC} -Wno-unused-result

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<
