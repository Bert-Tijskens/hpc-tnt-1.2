import numpy as np
import matplotlib.pyplot as plt
import math

bandwidth = 9.5 # GB/s
peak_performance = 11.2 # GFlops/s
I_sweet_spot = peak_performance/bandwidth


flops_per_atom = 2055.5
bytes_per_atom = 376.5
flops_per_byte = flops_per_atom/bytes_per_atom
print(flops_per_byte)
performance = 2.39

plt.plot([0,I_sweet_spot],[0,peak_performance],'b-'
        ,[I_sweet_spot,6],[peak_performance,peak_performance],'b-'
        ,[flops_per_byte],[performance],'ro'
        )
plt.xlabel('Computational Intensity [Flops/byte]')
plt.ylabel('Performance [GFlops/s]')
plt.show()