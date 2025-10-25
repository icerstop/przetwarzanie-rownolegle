# ---------- konfiguracja ----------
TARGET     := program
SRC        := main.cu            # lub tsp_heuristic.cu, jeœli tak nazwa³eœ plik

# docelowa architektura GPU; zmieñ wg w³asnej karty (np. sm_90, sm_89…)
GPU_ARCH   ?= sm_86

# Flagi host-kompilatora (gcc/clang) wywo³ywanego spod nvcc
HOSTFLAGS  := -O3 -fopenmp -finput-charset=UTF-8 -fexec-charset=UTF-8

# Kompletny zestaw flag nvcc
NVFLAGS    := -arch=$(GPU_ARCH) -std=c++17 \
              -Xcompiler '$(HOSTFLAGS)'     \
              -lineinfo -Wno-deprecated-gpu-targets

# ---------- regu³y g³ówne ----------
all: $(TARGET)

$(TARGET): $(SRC)
	nvcc $(NVFLAGS) -o $@ $<

run: $(TARGET)
	./$(TARGET) cities.csv 100 150 250

# Proste profilowanie GPU (wymaga Nsight Compute)
profile: $(TARGET)
	ncu --set full --target-processes all ./$(TARGET) cities.csv 100

clean:
	$(RM) $(TARGET) *.o *.nvvp

.PHONY: all clean run profile
