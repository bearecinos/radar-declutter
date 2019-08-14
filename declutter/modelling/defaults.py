'''Stores the default choices for wave, backscatter and directivity models.
This has to be in a separate module to "parameters.py" as that module contains
the variable which determines wavelength.
Were this class in that module, a cyclic import would occur when importing the
module. It would import "waves" to set the default wave, which would itself need
"parameters" to already be imported to define the wave.
There is likely a neater solution to this by restructuring where different
pieces of information are stored.

Having a class with methods to return the set option allows modules to
use None as the default and then load the chosen one from a central location.'''

import waves
import backscatter
import directivity

### should this also have an env instance to get time other parameters
### such as timestep or is additional nesting excessive?
class Defaults:
    wave = waves.GaussianDot()
    # static method treats function as value rather than
    # binding to instances of the class (can't return the function alone)
    intensity = staticmethod(backscatter.rayModel)
    directional = staticmethod(directivity.constant)
    def getWave(self):
        return self.wave
    def getIntensity(self):
        return self.intensity
    def getDirectivity(self):
        return self.directional

# like env in parameters, should only have 1 instance
default = Defaults()

# No wave to store a function so have to set each time module loaded
# If could store functions - could store a function from a module not
# normally imported and would get an error when function was unknown when read.
def setWave(wave):
    '''Excepts an instance of a wave class'''
    default.wave = wave
def setIntensity(intensity):
    default.intensity = intensity
def setDirectivity(directional):
    default.directional = directional


