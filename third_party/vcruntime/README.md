# Microsoft Visual C++ runtime 14.44.35211 (x64)

ONNX Runtime's Windows library (third_party/onnxruntime) needs Microsoft's
Visual C++ 2015-2022 runtime, which Windows does not have until something
installs it. These are its DLLs, unchanged, so the Windows program unpacks
them beside onnxruntime.dll and needs nothing installed (Microsoft allows
these files to be shipped with a program, "app-local").

Taken from Microsoft's vc_redist.x64.exe (https://aka.ms/vs/17/release/vc_redist.x64.exe,
SHA-256 cc0ff0eb1dc3f5188ae6300faef32bf5beeba4bdd6e8e445a9184072096b713b on
2026-11-21), cab a12, the *_amd64 files renamed.

| File | SHA-256 |
|---|---|
| msvcp140.dll | 0f885b509a685d2bbfa652fed26b5fb31d88fbdab0a978c641d1c7b8aa460aa9 |
| msvcp140_1.dll | bfad5aef4c63a669e3c140655cdfdf395b6c979b400a447bd5dcb65ed8826c3d |
| vcruntime140.dll | d5e4d9a3e835fa679450145d6a7d94e36573a509317111904d9b3712c30d9066 |
| vcruntime140_1.dll | 1f2d41c4aa5db0bc33ebf7b66d72943a817d7ce6cbe880502a9403823633093f |
