# ONNX Runtime's C API headers

`onnxruntime_c_api.h` and the two headers it includes, from Microsoft's official ONNX Runtime 1.31.1 release
(`onnxruntime-linux-x64-1.31.1.tgz`, `include/`; MIT License, `LICENSE`). Only the headers: `src/ft8_nn.cpp` loads
ONNX Runtime's library at run time (`dlopen` / `LoadLibrary`) when a model is given (`jt9 --nn`), so
nothing needs it otherwise; `jt9` looks for `libonnxruntime.so.1` (`onnxruntime.dll` on Windows)
beside itself or on the library path.
