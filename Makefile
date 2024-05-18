CXX:=g++

CXX_FLAGS_SHARED:=-Wall -Wextra 
CXX_FLAGS_RELEASE:=-O3
CXX_FLAGS_DEBUG:=-g

# ROOT dependency
ifeq (1, $(shell echo ${RHM_HAS_ROOT}))
	CXX_FLAGS_SHARED:=$(CXX_FLAGS_SHARED) `root-config --cflags --glibs`
else ifeq (, $(shell echo ${RHM_HAS_ROOT}))
	ifneq (, $(shell echo ${ROOTSYS}))
		CXX_FLAGS_SHARED:=$(CXX_FLAGS_SHARED) `root-config --cflags --glibs`
		RHM_HAS_ROOT:=1
	else
		RHM_HAS_ROOT:=0
	endif
endif

CXX_FLAGS_SHARED:=$(CXX_FLAGS_SHARED) -DRHM_HAS_ROOT=$(RHM_HAS_ROOT)

# Whether inputfile is expected
RHM_REQUIRES_INPUTFILE:=1

CXX_FLAGS_SHARED:=$(CXX_FLAGS_SHARED) -DRHM_REQUIRES_INPUTFILE=$(RHM_REQUIRES_INPUTFILE)

# Compile

HEADERS:=$(wildcard include/*.h)

APP_NAME = ${app}
SOURCES:=$(wildcard src/*.cxx)

OBJECTS:=$(SOURCES:.cxx=.o)
DEPENDENCIES:=$(OBJECTS:.o=.d)
PROGRAMS:=$(wildcard project/Cooling/*.cxx) $(wildcard project/M-R_diagram/*.cxx) $(wildcard project/Utility/*.cxx)

all : $(PROGRAMS)

$(PROGRAMS) : $(OBJECTS)
	@mkdir -p $(dir bin/${@:.cxx=.out})
	$(CXX) $@ $^ -o bin/${@:.cxx=.out} $(CXX_FLAGS_RELEASE) $(CXX_FLAGS_SHARED)

release : $(OBJECTS)
	@mkdir -p $(dir bin/$(APP_NAME))
	$(CXX) $(APP_NAME).cxx $^ -o bin/$(APP_NAME).out $(CXX_FLAGS_RELEASE) $(CXX_FLAGS_SHARED)

debug : $(OBJECTS)
	@mkdir -p $(dir bin/$(APP_NAME))
	$(CXX) $(APP_NAME).cxx $^ -o bin/$(APP_NAME).out $(CXX_FLAGS_DEBUG) $(CXX_FLAGS_SHARED)

%.o : %.cxx 
	$(CXX) -MMD -c $< -o $@

-include $(DEPENDENCIES)

clean:
	rm -f src/*.o src/*.d
	rm -r bin/*