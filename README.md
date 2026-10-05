# pd2af

[![License](https://img.shields.io/github/license/adrienrougny/pd2af)](https://github.com/adrienrougny/pd2af/blob/main/COPYING)

A library and CLI to transform process-description (PD) maps into activity-flow (AF) maps.

## Features

- **CellDesigner and SBGN-PD input**: a CellDesigner map gives a CellDesigner map, an SBGN-PD map gives an SBGN-AF map
- **Three transformation modes**: `normal`, `no-complex` and `keep-reactions`, with the `--keep-species` and `--drop-compartments` options
- **Three layout modes**: Graphviz `dot` auto-layout, reuse of the original positions, or an overlay on the original layout
- **Annotations and notes carried over**: from the input elements to the elements they produce

## Installation

```bash
pip install git+https://github.com/adrienrougny/pd2af
```

The `dot` layout mode requires Graphviz's `dot` binary on your `PATH`.

## Quick example

```python
import momapy.io.core
import pd2af

pd_map = momapy.io.core.read("my_map.xml").obj
af_map = pd2af.transform(pd_map).obj
momapy.io.core.write(af_map, "my_map_af.xml", writer="celldesigner")
```

```bash
pd2af transform my_map.xml -o my_map_af.xml
pd2af list-modes
```
