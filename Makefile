CXX:=g++

CXX_FLAGS_RELEASE:=-O3
CXX_FLAGS_DEBUG:=-Wall -Wextra -g
CXX_FLAGS_EXTRALIBS:=

# ROOT dependency
ifeq (1, $(shell echo ${RHM_HAS_ROOT}))
	CXX_FLAGS_EXTRALIBS:=$(CXX_FLAGS_EXTRALIBS) `root-config --cflags --glibs`
else ifeq (, $(shell echo ${RHM_HAS_ROOT}))
	ifneq (, $(shell echo ${ROOTSYS}))
		CXX_FLAGS_EXTRALIBS:=$(CXX_FLAGS_EXTRALIBS) `root-config --cflags --glibs`
		RHM_HAS_ROOT:=1
	else
		RHM_HAS_ROOT:=0
	endif
endif

# Whether inputfile is expected
RHM_REQUIRES_INPUTFILE:=1

HEADERS:=$(wildcard include/*.h)

APP_NAME = ${app}
SOURCES:=$(wildcard src/*.cxx)

OBJECTS:=$(SOURCES:.cxx=.o)
DEPENDENCIES:=$(OBJECTS:.o=.d)

all : release

release : $(OBJECTS)
	@mkdir -p $(dir bin/$(APP_NAME))
	$(CXX) $(APP_NAME).cxx $^ -o bin/$(APP_NAME).out $(CXX_FLAGS_RELEASE) $(CXX_FLAGS_EXTRALIBS) -DRHM_HAS_ROOT=$(RHM_HAS_ROOT) -DRHM_REQUIRES_INPUTFILE=$(RHM_REQUIRES_INPUTFILE)

debug : $(OBJECTS)
	@mkdir -p $(dir bin/$(APP_NAME))
	$(CXX) $(APP_NAME).cxx $^ -o bin/$(APP_NAME).out $(CXX_FLAGS_DEBUG) $(CXX_FLAGS_EXTRALIBS) -DRHM_HAS_ROOT=$(RHM_HAS_ROOT) -DRHM_REQUIRES_INPUTFILE=$(RHM_REQUIRES_INPUTFILE)

%.o : %.cxx 
	$(CXX) -MMD -c $< -o $@

-include $(DEPENDENCIES)

clean:
	rm -f src/*.o src/*.d
	rm -f bin/$(APP_NAME).out