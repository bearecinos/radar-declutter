# Get frequency/time step of radargram from models
from declutter.modelling import parameters
import numpy as np


# Template class
# offset - shifts wave earlier relative to response time on radargram
# only needed to allow IDL approach to be used too.
class Wave:
    offset = 0.0

    def amplitude(self, t):
        return 0.0


class GaussianDot(Wave):  # peak frequency equal to _freq
    freq = parameters.freq
    delt = 2.0/(3.0*freq)
    offset = 0.0

    def amplitude(self, t):
        return (-np.exp(0.5)*2*np.pi*self.freq*(t-self.delt) *
                np.exp(-2*np.pi**2*self.freq**2*(t-self.delt)**2))


# looks poor for one dataset but good for other?
class Ricker(Wave):  # peak frequency equal to _freq
    freq = parameters.freq
    delt = 2.0/(3.0*freq)
    offset = 0.0

    def amplitude(self, t):
        t = t - self.delt
        return (1-2*np.pi*np.pi*t*t*self.freq*self.freq)*np.exp(
            -np.pi*np.pi*self.freq*self.freq*t*t)


class Constant(Wave):  # equivalent to all frequencies, equal amplitude
    dt = parameters.env.getDt()
    delt = 0.0
    offset = 0.0

    def amplitude(self, t):
        return (t < self.dt) & (t >= 0)


class IDLWave(Wave):
    # exp( - (findgen(echo_size)-float(i))^2/50. )
    # timestep was 240us (2.4e-4 seconds echo length) / 2048 granularity
    # = 1.172e-7 = 117ns - much longer than our sampling rate
    # scaling = 1.172e-7
    dt = parameters.env.getDt()
    offset = parameters.env.getMaxTime() / 2.0
    delt = offset

    def __init__(self, c=50.0):
        self.c = float(c)

    def amplitude(self, t):
        t = (t-self.delt)/self.dt
        return np.exp(-t**2/self.c)


class RC(Wave):
    offset = 0.0

    def __init__(self, c=5.0):
        self.c = float(c)

    def amplitude(self, t):
        result = np.full_like(t, 0.0)  # get around vectorizing
        result[t > 0] = np.exp(-t[t > 0] / (self.c*1e-6))
        return result


class CosExp(Wave):
    freq = parameters.freq
    delt = 1.1/freq
    offset = 0.0

    def amplitude(self, t):
        p = 2*np.pi*self.freq*(t-self.delt)
        return np.cos(1.5*p)*np.exp(-p*p/16.0)


class Sym(Wave):
    freq = parameters.freq
    delt = 0.5/freq
    offset = 0.0

    def amplitude(self, t):
        p = 2*np.pi*self.freq*(t-self.delt)
        return (160*p-56*p*p*p+3*p**5)*np.exp(-p*p/4)/(1.343*64.0)
