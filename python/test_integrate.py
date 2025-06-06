# Run: `uv run test_integrate.py`
#
# /// script
# dependencies = [
#    'numpy',
#    'scipy',
#  ]
# ///

import numpy as np
from scipy.integrate import ode
import time
    
mu = 1000
y0 = np.r_[0.5, 0.5]
T = 500
tt = np.linspace(T/200,T, 100)

def van_der_pol(t, y):
    return np.r_[y[1], mu*(1.0-y[0]**2)*y[1]-y[0]]

c2 = time.process_time()
r2 = ode(van_der_pol).set_integrator('lsoda')
r2.set_initial_value(y0)
sol2 = np.array([r2.integrate(T) for T in tt])

for t, y in zip(tt, sol2[:,0]):
    print(t, y)
