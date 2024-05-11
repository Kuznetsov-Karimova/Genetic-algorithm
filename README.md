# Genetic-algorithm

University laboratory work


### How to setup and run
1. Ensure cmake and clang are installed properly
```bash
clang++ --version
```
```bash
cmake --version
```
```bash
make --version
```
2. Configure project:
```bash
cmake -S . -B build -G"MinGW Makefiles" -D"CMAKE_MAKE_PROGRAM:PATH=C:\Program Files (x86)\GnuWin32\bin\make"
```
Specify last flag with your make path. 


3. Build:
```bash
cmake --build build
```
