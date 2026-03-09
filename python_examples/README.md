# Python Examples

Python equivalents of the C++ examples under `examples/`.

## Build bindings first

```bash
cmake -S . -B build -DBUILD_PYTHON_BINDINGS=ON -DPython_EXECUTABLE=$(which python3)
cmake --build build --target sai_model_py -j8
```

## Run an example

```bash
python3 python_examples/01-create_model_from_file/main.py
```
