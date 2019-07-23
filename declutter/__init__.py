import warnings
# ignore warnings when comparing arrays containing NaNs.
# Always get False from ar > x where NaNs present
# which is desired behaviour.
warnings.filterwarnings('ignore')

def loadArcpy():
    import declutter.Generate
    import declutter.slopeAspect

__all__ = ["changeSize","models","path","stateless","test","loadArcpy"]    
