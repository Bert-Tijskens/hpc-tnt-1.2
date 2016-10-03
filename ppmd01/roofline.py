import numpy as np
import matplotlib.pyplot as plt
import math

yscale =1e9
peak_performance = 11.2e9/yscale
bandwidth        = 10.6e9/yscale
CI0 = peak_performance/bandwidth
flops = 14
bytes = 24
CI = flops/bytes
print('CI=',CI)
performance_random = bandwidth*CI/18.9

plt.plot( [0,CI0], [peak_performance,peak_performance], '--b'
        , [0,CI0], [0               ,peak_performance], '-b'
        , [CI0,5], [peak_performance,peak_performance], '-b'
        , [CI]   , [bandwidth*CI]                     , 'or'
        , [CI*24/28], [performance_random]            , 'or'
        )
plt.xlabel('Computational intensity [flops/byte]')
plt.ylabel('performance [Gflops/s]')
plt.title('Roofline')
plt.show()
