CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic -Iinclude -Ithird_party/eigen3

BUILD_DIR := build
LIB := $(BUILD_DIR)/libicp_registration.a
SRCS := src/kdtree.cpp src/rigid_transform.cpp src/icp.cpp
OBJS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

.PHONY: all test example clean

all: $(LIB) $(BUILD_DIR)/icp_example $(BUILD_DIR)/icp_tests

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: src/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIB): $(OBJS)
	ar rcs $@ $^

$(BUILD_DIR)/icp_tests: tests/test_icp.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $< $(LIB) -o $@

$(BUILD_DIR)/icp_example: examples/example_scan.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $< $(LIB) -o $@

test: $(BUILD_DIR)/icp_tests
	./$(BUILD_DIR)/icp_tests

example: $(BUILD_DIR)/icp_example
	./$(BUILD_DIR)/icp_example

clean:
	rm -rf $(BUILD_DIR)
