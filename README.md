# Area- and Topology-Preserving Polygon Simplification

## Overview
This project implements an **Area-Preserving Segment Collapse (APSC)** algorithm in **C++17** to simplify polygon boundaries while preserving:

- Exact **area** (within floating-point tolerance)
- **Topology** (no self-intersections, no ring crossings, same number of rings)
- Minimal **areal displacement**

The implementation follows the method described in *Kronenfeld et al. (2020)* and satisfies the requirements of **CSD2183 Data Structures Project 2**.

---

## Features
- Supports polygons with:
  - One exterior ring
  - Multiple interior rings (holes)
- Preserves:
  - Area of each ring
  - Overall polygon topology
- Minimizes areal displacement using APSC
- Efficient local updates after each collapse
- Output strictly follows assignment format

---

## Build Instructions

### Requirements
- g++ (C++17 or later)
- make

### Build
Run in project root:
```
make
```

This generates:
```
./simplify
```

---

## Usage

```
./simplify <input_file.csv> <target_vertices>
```

### Example
```
./simplify tests/input_rectangle_with_two_holes.csv 7
```

---

## Input Format

CSV format:
```
ring_id,vertex_id,x,y
```

Notes:
- ring_id = 0 → exterior ring (counterclockwise)
- ring_id > 0 → interior rings (holes, clockwise)
- Input is assumed valid and non-self-intersecting

---

## Output Format

Program prints to standard output:
```
ring_id,vertex_id,x,y
...
Total signed area in input: <value>
Total signed area in output: <value>
Total areal displacement: <value>
```

Requirements:
- Same ring structure as input
- Vertex IDs start from 0 and are contiguous
- Values printed in scientific notation
- No extra stdout output (use stderr for debug)

---

## Algorithm Summary

The APSC algorithm works as follows:

1. Consider sequences of four vertices: A → B → C → D  
2. Replace with A → E → D  
3. Compute E such that:
   - Area is preserved exactly
   - Areal displacement is minimized  
4. Select the collapse with minimum displacement  
5. Validate topology (no intersections)  
6. Update locally  
7. Repeat until target vertex count is reached  

---

## Data Structures

- Ring-based polygon representation  
- Local candidate evaluation (A, B, C, D sequences)  
- Incremental updates after each collapse  
- Intersection checks for topology preservation  

---

## Tests

Test files are located in:
```
tests/
```

Run all tests:
```
make test
```

---

## Test Results


All outputs:
- Near expected results  
- Preserve area within tolerance  
- Maintain topology  

---

## Performance

- Uses local updates instead of full recomputation  
- Suitable for large inputs (≥100k vertices)  
- Approximate complexity: O(n log n) depending on implementation  

---

## Limitations

- Minor floating-point precision errors possible  
- Very constrained inputs may limit further simplification  

---

## References

Kronenfeld, B. J. et al. (2020)  
*Simplification of polylines by segment collapse: minimizing areal displacement while preserving area*

---

## Notes

- Focus on correctness, efficiency, and clean design
