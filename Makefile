CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -pedantic

SRC = src/main.cpp src/geometry.cpp src/polygon.cpp src/simplify.cpp
OUT = simplify
TEST_DIR = tests
OUTPUT_DIR = $(TEST_DIR)/output
INPUTS = $(wildcard $(TEST_DIR)/input_*.csv)

.PHONY: all build test clean clean-test

all: $(OUT) test

build: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT)

test: $(OUT)
	@mkdir -p "$(OUTPUT_DIR)"; \
	total=0; pass_count=0; fail_count=0; \
	fail=0; \
	for input in $(sort $(INPUTS)); do \
		total=$$((total+1)); \
		base=$$(basename "$$input"); \
		name=$${base#input_}; \
		name=$${name%.csv}; \
		expected="$(TEST_DIR)/output_$$name.txt"; \
		actual="$(OUTPUT_DIR)/$$name.txt"; \
		case "$$name" in \
			original_*) target=99 ;; \
			rectangle_with_two_holes) target=7 ;; \
			cushion_with_hexagonal_hole) target=13 ;; \
			blob_with_two_holes) target=17 ;; \
			wavy_with_three_holes) target=21 ;; \
			lake_with_two_islands) target=17 ;; \
			*) target=99 ;; \
		esac; \
		printf "RUN %s\n" "$$name"; \
		./$(OUT) "$$input" "$$target" > "$$actual"; rc=$$?; \
		if [ $$rc -ne 0 ]; then \
			printf "FAIL %s (program exit %s)\n" "$$name" "$$rc"; \
			fail_count=$$((fail_count+1)); \
			fail=1; \
			continue; \
		fi; \
		if [ ! -f "$$expected" ]; then \
			printf "FAIL %s (missing expected output)\n" "$$name"; \
			fail_count=$$((fail_count+1)); \
			fail=1; \
		else \
			if diff -q "$$expected" "$$actual" >/dev/null; then \
				printf "PASS %s\n" "$$name"; \
				pass_count=$$((pass_count+1)); \
			else \
				printf "FAIL %s\n" "$$name"; \
				fail_count=$$((fail_count+1)); \
				printf "%s\n" "----- DIFF ($$name) BEGIN -----"; \
				diff -u "$$expected" "$$actual" || true; \
				printf "%s\n" "----- DIFF ($$name) END -----"; \
				fail=1; \
			fi; \
		fi; \
	done; \
	printf "\nTOTAL %s | PASS %s | FAIL %s\n" "$$total" "$$pass_count" "$$fail_count"; \
	exit $$fail

clean-test:
	rm -rf "$(OUTPUT_DIR)"

clean: clean-test
	rm -f $(OUT)
