import numpy as np
import matplotlib.pyplot as plt

r = np.linspace(.9,5,411)
ri2 = 1/(r*r)
ri6 = ri2*ri2*ri2
lj   = 4*ri6*(ri6-1)
ljff = 48*ri2*ri6*(ri6-0.5);
plt.plot(r,lj,r,ljff,[0.5,5],[0,0])
plt.axis([0.5,5,-2,1])
plt.legend(['lj(r)','lj_force_factor(r)'])
plt.show()

#         rr = 1.0d0/r2
#         rr6 = rr*rr*rr;
#         lj_pot2 = 4.0d0*rr6*(rr6-1.0d0);
