COMPILER := g++
CPPFLAGS := -O2 -std=c++17 -Wall -Wno-unused-variable
OBJECTS  := AMPC.cpp
DFLAGS = 

default_target: AMPC

AMPC:
	@$(COMPILER) $(CPPFLAGS) $(DFLAGS) $(OBJECTS) -o AMPC.exe

test:
	@$(COMPILER) $(CPPFLAGS) $(DFLAGS) test.cpp -o test.exe
