from numpy import *
from general_fit import *
from matplotlib.pyplot import *
from uncertainties import ufloat
from uncertainties.umath import *

voltages_error = 1.0
time_error = 2.0
def model_cube(parameters, data):
    a = parameters[0]
    b = parameters[1]
    c = parameters[2]
    d = parameters[3]
    
    return a*data**3 + b*data**2 + c*data + d

def residuals_cube(parameters, x, y):
    return y - model_cube(parameters, x)

figure(figsize = (11, 5))
suptitle("Scintillator 1")
### Scintillator 1
s1_run1_counts = array([244,402,644,777,936,1064,1340,1692,3439,8777])+0.0
s1_run1_times = ones_like(s1_run1_counts)*180.0
s1_run1_voltages = -array([1900,1920,1940,1960,1980,2000,2020,2040,2060,2080])+0.0

s1_run1_counts_per_second = s1_run1_counts/s1_run1_times
s1_run1_counts_per_second_error = s1_run1_counts * time_error / s1_run1_times**2

guess = polyfit(s1_run1_voltages, s1_run1_counts_per_second, 3).tolist()
odr_parameter_ideal, odr_parameter_error, odr_parameter_p_value, odr_SSR, odr_covariance = general_fit(s1_run1_voltages, s1_run1_counts_per_second, model_cube, guess, x_err = ones_like(s1_run1_voltages)*voltages_error, y_err = s1_run1_counts_per_second_error, return_covariance = True)

s1_r1_ssr = odr_SSR

subplot(1, 2, 1)
x = linspace(s1_run1_voltages[0], s1_run1_voltages[-1], 1000)
f = lambda z : model_cube(odr_parameter_ideal, z)
errorbar(s1_run1_voltages, s1_run1_counts_per_second, xerr = ones_like(s1_run1_voltages)*voltages_error, yerr = s1_run1_counts_per_second_error, color = 'r', marker = 'x', linestyle = '', label = 'Counts per second, first run')
plot(x, f(x), 'k--', label = 'Fitted cubic polynomial, first run')

a = odr_parameter_ideal[0]
sigma_a = odr_parameter_error[0]
b = odr_parameter_ideal[1]
sigma_b = odr_parameter_error[1]

working_point_s1_r1 = -b/(3*a)
working_point_s1_r1_error = sqrt(sigma_b**2/(9.0*a**2) + b**2*sigma_a**2/(9.0*a**4) - 2.0 * b/(9.0*a**3) * odr_covariance[0, 1] )

print "The working point of scintillator 1 for the first run is ", working_point_s1_r1, "+-", working_point_s1_r1_error, "V"
print "The p-value is: ", odr_parameter_p_value
s1_r1_p_value = odr_parameter_p_value

s1_r1_parameters = odr_parameter_ideal

s1_run2_counts = array([256,265,317,321,344,391,384,494,579,670]) + 0.0
s1_run2_times = ones_like(s1_run2_counts)*60.0
s1_run2_voltages = -array([1960,1970,1980,1990,2000,2010,2020,2030,2040,2050])+0.0

s1_run2_counts_per_second = s1_run2_counts / s1_run2_times
s1_run2_counts_per_second_error = s1_run2_counts * time_error / s1_run2_times**2

guess = polyfit(s1_run2_voltages, s1_run2_counts_per_second, 3).tolist()
odr_parameter_ideal, odr_parameter_error, odr_parameter_p_value, odr_SSR, odr_covariance = general_fit(s1_run2_voltages, s1_run2_counts_per_second, model_cube, guess, x_err = ones_like(s1_run2_voltages)*voltages_error, y_err = s1_run2_counts_per_second_error, return_covariance = True)

s1_r2_ssr = odr_SSR

subplot(1, 2, 1)
x = linspace(s1_run2_voltages[0], s1_run2_voltages[-1], 1000)
f = lambda z : model_cube(odr_parameter_ideal, z)
errorbar(s1_run2_voltages, s1_run2_counts_per_second, xerr = ones_like(s1_run2_voltages)*voltages_error, yerr = s1_run2_counts_per_second_error, color = 'b', marker = '.', linestyle = '', label = 'Counts per second, second run')
plot(x, f(x), 'k:', label = 'Fitted cubic polynomial, second run')
title('Plateau manifesting for the counts per second')
xlabel('Voltage in V')
ylabel('Counts per second in Hz')
legend()

a = odr_parameter_ideal[0]
sigma_a = odr_parameter_error[0]
b = odr_parameter_ideal[1]
sigma_b = odr_parameter_error[1]

working_point_s1_r2 = -b/(3*a)
working_point_s1_r2_error = sqrt(sigma_b**2/(9.0*a**2) + b**2*sigma_a**2/(9.0*a**4) - 2.0 * b/(9.0*a**3) * odr_covariance[0, 1] )

print "The working point of scintillator 1 for the second run is ", working_point_s1_r2, "+-", working_point_s1_r2_error, "V"
print "The p-value is: ", odr_parameter_p_value
s1_r2_p_value = odr_parameter_p_value

print "The mean for scintillator 1 is then {:4.0f}V".format((working_point_s1_r1 + working_point_s1_r2)/2.0)

s1_r2_parameters = odr_parameter_ideal

subplot(1, 2, 2)
plot(s1_run1_voltages, zeros_like(s1_run1_voltages), 'k')
plot(s1_run1_voltages, residuals_cube(s1_r1_parameters, s1_run1_voltages, s1_run1_counts_per_second), 'rx', linestyle = '', label = "First run. SSR = {}".format(s1_r1_ssr) + "\np-value = {}".format(s1_r1_p_value))
plot(s1_run2_voltages, residuals_cube(s1_r2_parameters, s1_run2_voltages, s1_run2_counts_per_second), 'b.', linestyle = '', label = "Second run. SSR = {}".format(s1_r2_ssr) + "\np-value = {}".format(s1_r2_p_value))
title('Residuals of fits')
ylabel('Error in the counts per second, in Hz')
xlabel('Voltage in V')
legend()
savefig('scintillator1.png')
clf()


figure(figsize = (11, 5))
suptitle("Scintillator 2")
### Scintillator 2
subplot(1, 2, 1)
s2_run1_counts = array([143,215,324,329,388,423,538,667])+0.0
s2_run1_times = ones_like(s2_run1_counts)*60.0
s2_run1_voltages = -array([1900,1920,1940,1960,1980,2000,2020,2040])+0.0

s2_run1_counts_per_second = s2_run1_counts/s2_run1_times
s2_run1_counts_per_second_error = s2_run1_counts * time_error / s2_run1_times**2

guess = polyfit(s2_run1_voltages, s2_run1_counts_per_second, 3).tolist()
odr_parameter_ideal, odr_parameter_error, odr_parameter_p_value, odr_SSR, odr_covariance = general_fit(s2_run1_voltages, s2_run1_counts_per_second, model_cube, guess, x_err = ones_like(s2_run1_voltages)*voltages_error, y_err = s2_run1_counts_per_second_error, return_covariance = True)

s2_r1_ssr = odr_SSR

x = linspace(s2_run1_voltages[0], s2_run1_voltages[-1], 1000)
f = lambda z : model_cube(odr_parameter_ideal, z)
errorbar(s2_run1_voltages, s2_run1_counts_per_second, xerr = ones_like(s2_run1_voltages)*voltages_error, yerr = s2_run1_counts_per_second_error, color = 'r', marker = 'x', linestyle = '', label = 'Counts per second, first run')
plot(x, f(x), 'k--', label = 'Fitted cubic polynomial, first run')

a = odr_parameter_ideal[0]
sigma_a = odr_parameter_error[0]
b = odr_parameter_ideal[1]
sigma_b = odr_parameter_error[1]

working_point_s2_r1 = -b/(3*a)
working_point_s2_r1_error = sqrt(sigma_b**2/(9.0*a**2) + b**2*sigma_a**2/(9.0*a**4) - 2.0 * b/(9.0*a**3) * odr_covariance[0, 1] )

print "The working point of scintillator 2 for the first run is ", working_point_s2_r1, "+-", working_point_s2_r1_error, "V"
print "The p-value is: ", odr_parameter_p_value

s2_r1_p_value = odr_parameter_p_value

s2_r1_parameters = odr_parameter_ideal

s2_run2_counts = array([697,895,860,871,1047,1219]) + 0.0
s2_run2_times = array([180,210,180,180,180,180])+0.0
s2_run2_voltages = -array([1920,1930,1940,1950,1960,1970])+0.0

s2_run2_counts_per_second = s2_run2_counts / s2_run2_times
s2_run2_counts_per_second_error = s2_run2_counts * time_error / s2_run2_times**2

guess = polyfit(s2_run2_voltages, s2_run2_counts_per_second, 3).tolist()
odr_parameter_ideal, odr_parameter_error, odr_parameter_p_value, odr_SSR, odr_covariance = general_fit(s2_run2_voltages, s2_run2_counts_per_second, model_cube, guess, x_err = ones_like(s2_run2_voltages)*voltages_error, y_err = s2_run2_counts_per_second_error, return_covariance = True)

s2_r2_parameters = odr_parameter_ideal
s2_r2_ssr = odr_SSR

x = linspace(s2_run2_voltages[0], s2_run2_voltages[-1], 1000)
f = lambda z : model_cube(odr_parameter_ideal, z)
errorbar(s2_run2_voltages, s2_run2_counts_per_second, xerr = ones_like(s2_run2_voltages)*voltages_error, yerr = s2_run2_counts_per_second_error, color = 'b', marker = '.', linestyle = '', label = 'Counts per second, second run')
plot(x, f(x), 'k:', label = 'Fitted cubic polynomial, second run')
title('Plateau manifesting for the counts per second')
xlabel('Voltage in V')
ylabel('Counts per second in Hz')
legend()

a = odr_parameter_ideal[0]
sigma_a = odr_parameter_error[0]
b = odr_parameter_ideal[1]
sigma_b = odr_parameter_error[1]

working_point_s2_r2 = -b/(3*a)
working_point_s2_r2_error = sqrt(sigma_b**2/(9.0*a**2) + b**2*sigma_a**2/(9.0*a**4) - 2.0 * b/(9.0*a**3) * odr_covariance[0, 1] )

print "The working point of scintillator 2 for the first run is ", working_point_s2_r2, "+-", working_point_s2_r2_error, "V"
print "The p-value is: ", odr_parameter_p_value
s2_r2_p_value = odr_parameter_p_value

print "The mean for scintillator 2 is then {:4.0f}V".format((working_point_s2_r1 + working_point_s2_r2)/2.0)

subplot(1, 2, 2)
plot(s2_run1_voltages, zeros_like(s2_run1_voltages), 'k')
plot(s2_run1_voltages, residuals_cube(s2_r1_parameters, s2_run1_voltages, s2_run1_counts_per_second), 'rx', linestyle = '', label = "First run. SSR = {}".format(s2_r1_ssr) + "\np-value = {}".format(s2_r1_p_value))
plot(s2_run2_voltages, residuals_cube(s2_r2_parameters, s2_run2_voltages, s2_run2_counts_per_second), 'b.', linestyle = '', label = "Second run. SSR = {}".format(s2_r2_ssr) + "\np-value = {}".format(s2_r2_p_value))
title('Residuals of fits')
ylabel('Error in the counts per second, in Hz')
xlabel('Voltage in V')
legend()
savefig('scintillator2.png')
clf()




