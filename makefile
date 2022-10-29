# Makefile for C++

COMPILER = g++
VERSION = -std=c++14
ABSOLUTE_PROJECT_DIR = $(shell pwd)
ENTRY_NAME = Application

ABSOLUTE_DEPENDENCY_PATH_1 = $(ABSOLUTE_PROJECT_DIR)/include/Math
# ABSOLUTE_DEPENDENCY_PATH_2 = $(ABSOLUTE_PROJECT_DIR)/libs/Utility

ASM_FLAGS = $(VERSION) -S -fverbose-asm -g
BIN_FLAGS = $(VERSION) -Wall -Wextra -g
OBJDUMP_FLAGS = -S --disassemble
DEPENDENCY_FLAGS = -I$(ABSOLUTE_DEPENDENCY_PATH_1)
# DEPENDENCY_FLAGS = -I$(ABSOLUTE_DEPENDENCY_PATH_1) -I$(ABSOLUTE_DEPENDENCY_PATH_2)

# $(COMPILER) $(ASM_FLAGS) -o build/$(ENTRY_NAME).asm $(ENTRY_NAME).cpp
# $(info $$ASM_FLAGS is [${ASM_FLAGS}])
# include/Math/Algebra/DualNumber.cpp
# %.o: %.cpp
# 	$(COMPILER) $(BIN_FLAGS) $(DEPENDENCY_FLAGS) -o $@

# $(ENTRY_NAME): $(ENTRY_NAME).cpp
# 	mkdir -p build
# 	$(shell make clean)
# 	$(COMPILER) $(BIN_FLAGS) $(DEPENDENCY_FLAGS) -o build/$(ENTRY_NAME) $(ENTRY_NAME).cpp include/Math/Algebra/DualNumber.cpp
# 	$(COMPILER) $(ASM_FLAGS) $(DEPENDENCY_FLAGS) -o build/$(ENTRY_NAME).asm $(ENTRY_NAME).cpp
# 	objdump $(OBJDUMP_FLAGS) build/$(ENTRY_NAME) > build/$(ENTRY_NAME).dump
# 	@echo "Build successfully.\n"


$(ENTRY_NAME): $(ENTRY_NAME).cpp
	mkdir -p build
	$(shell make clean)
	$(COMPILER) $(BIN_FLAGS) $(DEPENDENCY_FLAGS) -o build/$(ENTRY_NAME) $(ENTRY_NAME).cpp
	$(COMPILER) $(ASM_FLAGS) $(DEPENDENCY_FLAGS) -o build/$(ENTRY_NAME).asm $(ENTRY_NAME).cpp
	objdump $(OBJDUMP_FLAGS) build/$(ENTRY_NAME) > build/$(ENTRY_NAME).dump
	@echo "Build successfully.\n"


clean:
	rm -f build/$(ENTRY_NAME)
	rm -f build/*.asm
	rm -f *.s
	rm -f build/$(ENTRY_NAME).dump