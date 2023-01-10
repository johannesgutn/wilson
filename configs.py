import numpy as np
from scipy.integrate import odeint
import math as mt


# Values for parameters
Nc=3
z=0.5
qhatmix=1.5 #GeV^2/fm
qhat=qhatmix*25.77 #fm^(-3)
theta=0.5
CF=(Nc**2-1)/(2*Nc)
t1=0.3
t2=1
tau = t2-t1
a = qhat*theta**2/(2*Nc)
t_max=t2+3


eGev = 1000 # Energy in Gev
e = eGev * 5.076 # Conversion factor to fm

# Names for saving
filename = 'wilson_line'


n1=np.array([0,(1-z)*theta])
n2=np.array([0,-z*theta])

x=lambda t: (t-t1)*n1
y=lambda t:(t-t1)*n2
X=lambda t:(t-t2)*n1
Y=lambda t:(t-t2)*n2

r=[x,Y]
R=[y,X]