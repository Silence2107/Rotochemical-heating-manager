CXX:=g++

CXXFLAGS:=-Wall -Wextra 
all: CXXFLAGS +=-O3
release: CXXFLAGS +=-O3
debug: CXXFLAGS +=-g

# ROOT dependency
ifeq (1, $(shell echo ${RHM_HAS_ROOT}))
	CXXFLAGS += `root-config --cflags --glibs`
else ifeq (, $(shell echo ${RHM_HAS_ROOT}))
	ifneq (, $(shell echo ${ROOTSYS}))
		CXXFLAGS += `root-config --cflags --glibs`
		RHM_HAS_ROOT:=1
	else
		RHM_HAS_ROOT:=0
	endif
endif

CXXFLAGS += -DRHM_HAS_ROOT=$(RHM_HAS_ROOT)

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
	$(CXX) $@ $^ -o bin/${@:.cxx=.out} $(CXXFLAGS)

release : $(OBJECTS)
	@mkdir -p $(dir bin/$(APP_NAME))
	$(CXX) $(APP_NAME).cxx $^ -o bin/$(APP_NAME).out $(CXXFLAGS)

debug : $(OBJECTS)
	@mkdir -p $(dir bin/$(APP_NAME))
	$(CXX) $(APP_NAME).cxx $^ -o bin/$(APP_NAME).out $(CXXFLAGS)

%.o : %.cxx 
	$(CXX) -MMD -c $< -o $@

-include $(DEPENDENCIES)

clean:
	rm -f src/*.o src/*.d
	rm -r bin/*