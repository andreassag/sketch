CXX ?= g++
CXXFLAGS ?= -std=c++20 -O3 -Wall -Wextra -Wpedantic -Wconversion -Iinclude

SRCDIR = src
INCDIR = include
BINDIR = bin

all: $(BINDIR)/sketch

$(BINDIR):
	mkdir -p $(BINDIR)

$(BINDIR)/sketch: $(SRCDIR)/main.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BINDIR)/sketch_tests: tests/test_sketch.cpp | $(BINDIR)
	$(CXX) $(CXXFLAGS) $< -o $@

test: $(BINDIR)/sketch_tests
	./$(BINDIR)/sketch_tests

clean:
	rm -rf $(BINDIR) build *.skh

.PHONY: all test clean
