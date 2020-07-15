import numpy as np
from scipy.optimize import curve_fit
import random
import matplotlib.pyplot as plt


def quadratic(x, q, m, c):
    """
    Quatratic function
    Inputs:
    #	data		:x (float)
    #	Constant	:c (float)
    #	linear para :m (float)
    #	quad part 	:q (float)
    """
    y = []
    for i in range(len(x)):
        y.append(q*x[i]*x[i] + m*x[i] + c)
    return y


def chi2(y_obs, y_err, y_exp, n_par):
    chi2 = 0
    ndof = len(y_obs) - n_par - 1
    for i in range(len(y_exp)):
        chi2 += ((y_exp[i] - y_obs[i])/y_err[i])**2
    chi2 = chi2/ndof
    return chi2


noise = 10*np.asarray(random.sample(range(-100, 100), 200))

x_array = np.linspace(-100,100,200)

y_array = np.array(quadratic(x_array, 2, 1, 0)) + noise

y_err = np.ones_like(x_array)*1000

p_guess = [1,1,1]
p_bounds = [[0,0,0],[10,10,10]]

popt, pcov = curve_fit(quadratic, x_array, y_array, p0=p_guess, bounds=p_bounds)
print("The optimised fitted parameters are: ", popt)
print("The covariance matrix is: ", pcov)
for i in range(len(popt)):
    print("Error on parameter {}: {} is {}".format(i, popt[i], np.sqrt(pcov[i][i])))

print("The reduced chi2 is: ", chi2(y_array, y_err, quadratic(x_array, *popt), 3))

plt.errorbar(x_array, y_array, y_err, color='r', fmt='.')
plt.plot(x_array, quadratic(x_array, *popt), color='b')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
