# Models for `jt9 --nn`

`ft8nn-jt9.onnx` (3.5 MB): the neural FT8 soft demodulator for `jt9`. Use it as
`jt9 -8 --nn models/ft8nn-jt9.onnx` (optionally `--nn-level 2`), with ONNX Runtime's library
(`libonnxruntime.so.1`, `onnxruntime.dll` on Windows; its official release, 1.31.1 or newer) beside `jt9` or on the
library path: `jt9` loads it only for `--nn`. Each release carries both: `ft8nn-jt9.onnx`, and the library in
`onnxruntime-linux-x86_64.tgz`, `onnxruntime-linux-aarch64.tgz` or `onnxruntime-windows-x86_64.zip` (Microsoft's
1.31.1, unchanged, with its licence), to unpack beside `jt9`.

- Network: a transformer (4 layers, d 128, ~850k parameters); input `features` [B, 79, 48],
  output `llr` [B, 174] (log P(1)/P(0)), runs on the CPU through ONNX Runtime.
