## Installation
```bash
apt-get update && apt-get install -y liblapack-dev liblapacke-dev
```

```bash
pip install -e .
```

## Usage
```python
from eiqp import EIQP
import numpy as np

# Initialize solver
solver = EIQP()

# Setup your problem
Q = np.eye(6)
c = np.array([[-1.0], [-0.8], [-0.6], [-0.4], [-0.2], [0.0]])
A = -np.eye(6)
b = -np.ones(6)

# Solve
z, status = solver.solve(Q, c, A, b, epsilon=1e-8)

print("Status:", status)
print("Solution:", z)
# Status: True
# Solution: [9.99980090e-01 7.99999998e-01 5.99999999e-01 4.00000000e-01
#  2.00000001e-01 1.99100974e-05]
```