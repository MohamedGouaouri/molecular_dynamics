LIB        = -L. 
INCLUDE    = -I.
CFLAGS     = -O3
EXEC       = MD
BUILD_PATH = build
CXX        = g++

${EXEC}: MD.cpp
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} MD.cpp -o ${BUILD_PATH}/${EXEC}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<
