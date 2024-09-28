This is the code repository for our paper "[Optimizing T gates in Clifford+T circuit as $π/4$ rotations around Paulis](https://arxiv.org/abs/1903.12456)".

This code depends on [QuaEC](https://github.com/cgranade/python-quaec), which can be installed by cloning the QuaEC repo and running `pip install .` (directly doing `pip install QuaEC` doesn't seem to work in our environment). If Python complains about `from collections import Sequence`, change it to `from collections.abc import Sequence`

Other dependencies `numpy`, `gmpy2`, and `Cython` (for `pyximport`) can be installed with pip directly.

To run tests:
```
python3 -m optimize.benchmark <circuit file>
```
The circuit files are in the `benchmark` folder. Note that some circuits in the `Nam` subfolder are not supported because they contain `RZ` gates of arbitrary angles.
