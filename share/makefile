DSOURCES := ${wildcard src/*.d}
CSOURCES := src/beta.o

LDC := ldc
GDC := gdc
DMD := dmd

CHECK_LDC := ${shell command -v ${LDC} 2> /dev/null}
CHECK_GDC := ${shell command -v ${GDC} 2> /dev/null}
CHECK_DMD := ${shell command -v ${DMD} 2> /dev/null}


ifneq (${CHECK_LDC},)
	COMPILER := ${LDC}
	RELEASE_FLAGS := -Jviews -release -enable-inlining -O -w -oq
	DEBUG_FLAGS := -Jviews -d-debug -g -unittest -w
	OUTPUT_FLAG := -of
else
ifneq (${CHECK_GDC},)
	COMPILER := ${GDC}
	RELEASE_FLAGS := -Jviews -frelease -finline-functions -O3 -Werror -Wall
	DEBUG_FLAGS := -Jviews -fdebug -g -funittest -Werror -Wall
	OUTPUT_FLAG := -o
	STATIC_FLAGS := -lgsl -lgslcblas
else
	COMPILER := ${DMD}
	RELEASE_FLAGS := -Jviews -release -inline -O -noboundscheck
	DEBUG_FLAGS := -Jviews -debug -g -unittest -w
	OUTPUT_FLAG := -of
endif
endif

ifeq (${CHECK_LDC},)
ifeq (${CHECK_LDC},)
ifeq (${CHECK_DMD},)
${error No D compiler found at ${LDC} or ${DMD} or ${GDC}}
endif
endif
endif

CLEAN_OBJECTS := rm -f src/*.o bin/*.o *.o

bin/alfonso	: ${DSOURCES} ${CSOURCES} views/commit
	${COMPILER} ${RELEASE_FLAGS} ${DSOURCES} ${CSOURCES} ${STATIC_FLAGS} ${OUTPUT_FLAG}bin/alfonso
	${CLEAN_OBJECTS}

test	: ${DSOURCES} ${CSOURCES} views/commit
	${COMPILER} ${DEBUG_FLAGS} ${DSOURCES} ${CSOURCES} ${STATIC_FLAGS} ${OUTPUT_FLAG}unittest
	./unittest
	${CLEAN_OBJECTS} unittest

views/commit : ${DSOURCES} ${CSOURCES}
	mkdir -p views
	git rev-parse --short HEAD > views/commit

.PHONY : test clean

clean :
	${CLEAN_OBJECTS} bin/alfonso

