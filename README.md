# LSAL Accelerator — Smith-Waterman Local Sequence Alignment

Hardware/software co-design project implementing the **Smith-Waterman** local sequence alignment algorithm across multiple platforms and optimization levels. Developed as part of an embedded systems / hardware acceleration course.

## Overview

Smith-Waterman is a classic dynamic-programming algorithm for finding the best local alignment between two DNA sequences. Given a query `Q` and a database sequence `D`, it fills an `M×N` scoring matrix and then traces back through the matrix to recover the aligned subsequences.

Scoring parameters used throughout:
| Event    | Score |
|----------|-------|
| Match    | +2    |
| Mismatch | −1    |
| Gap      | −1    |

## Repository Structure

```
.
├── x86/                    # Implementations for x86 CPU
│   ├── lsal_u_x86.c        # Unoptimized baseline (flat linear index)
│   ├── lsal_o_x86.c        # Optimized (row-major nested loops)
│   └── lsal_omp_x86.c      # Parallel — OpenMP wavefront tiling
│
├── arm/                    # Implementations for ARM CPU
│   ├── lsal_u_arm.c        # Unoptimized baseline
│   ├── lsal_opt_arm.c      # Optimized (row-major nested loops)
│   └── lsal_par_arm.c      # Parallel — anti-diagonal wavefront
│
├── hls/                    # FPGA accelerator (Xilinx Vitis HLS)
│   ├── lsal.h              # Kernel function declaration
│   ├── lsal.cpp            # HLS kernel with AXI/pipeline pragmas
│   └── lsal_host.cpp       # OpenCL host code (runs on ARM of MPSoC)
│
└── Embbeded_report.pdf     # Full project report
```

## Implementations

### x86 / ARM — Three Optimization Levels

Each platform has three progressively optimized C implementations:

1. **Unoptimized (`_u`)** — Computes both `row` and `col` via integer division on a flat index. Simple but inefficient due to costly division operations.

2. **Optimized (`_o` / `_opt`)** — Uses explicit nested `row`/`col` loops and direct index arithmetic, improving cache locality and removing the division overhead.

3. **Parallel (`_omp` / `_par`)** — Exploits the anti-diagonal (wavefront) dependency structure of the DP matrix. Cells on the same anti-diagonal are independent and can be computed concurrently. The x86 version uses OpenMP with tile-level wavefront scheduling; the ARM version iterates anti-diagonals directly.

### FPGA Accelerator (Xilinx Vitis HLS + OpenCL)

The HLS kernel (`lsal.cpp`) is designed for a fixed query length `N = 32` and database length `M = 65536`. Key design decisions:

- **Anti-diagonal streaming** — The database is padded and streamed through a sliding window, so all `N` query columns are processed in parallel each clock cycle.
- **`#pragma HLS PIPELINE`** — The outer `Round` loop is pipelined so the FPGA issues one anti-diagonal per clock.
- **Complete array partitioning** — `q_buf`, score buffers, and direction buffers are fully partitioned, giving simultaneous access to all `N` elements.
- **AXI interfaces** — Query and database are accessed over AXI master (`m_axi`) ports; `max_idx` is returned over AXI-Lite (`s_axilite`).

The host code (`lsal_host.cpp`) runs on the ARM cores of the Zynq MPSoC (e.g., ZCU102) and drives the FPGA kernel via the OpenCL API. After the hardware run it re-executes the reference software implementation and compares results for verification.

## Building

### x86 (GCC)

```bash
# Optimized sequential
gcc -O2 -o lsal_o x86/lsal_o_x86.c

# Parallel (OpenMP)
gcc -O2 -fopenmp -o lsal_omp x86/lsal_omp_x86.c

# Run: <query_length> <database_length>
./lsal_o 128 1024
./lsal_omp 128 1024
```

Compile with `-DTEST=1` to enable traceback output and print the alignment instead of timing.

### ARM (cross-compile or native)

Same flags as x86. For cross-compilation, replace `gcc` with your ARM toolchain (e.g., `aarch64-linux-gnu-gcc`).

### FPGA (Xilinx Vitis)

Synthesis and place-and-route are performed with Vitis HLS targeting the board's part number. The output `.xclbin` binary is then passed to the host application:

```bash
./lsal_host <path/to/kernel.xclbin>
```

## Dependencies

| Component | Dependency |
|-----------|-----------|
| x86 parallel | GCC + OpenMP (`-fopenmp`) |
| ARM parallel | GCC |
| FPGA kernel  | Xilinx Vitis HLS, `ap_int.h` |
| FPGA host    | OpenCL (`CL/opencl.h`, `CL/cl_ext.h`), Xilinx runtime (XRT) |

## Report

See [`Embbeded_report.pdf`](Embbeded_report.pdf) for the full analysis, including design rationale, performance measurements across platforms, and discussion of the HLS optimization choices.
