import numpy as np
import matplotlib.pyplot as plt

mu = 0.0
sigma = 1.0
x = mu + sigma*np.random.normal(size=5000)

plt.figure()
plt.hist(x, bins=20, alpha=0.8)
plt.axvline(np.mean(x), ls='--', color='k', lw=2.0, label=r'$\mu$='+str(mu)+', $\sigma$ = '+str(sigma))
plt.axvline(0., ls='--', color='b', lw=2.0)
plt.legend()
plt.show()
