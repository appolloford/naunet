from .dust import Dust


class UniDust(Dust):

    vars = {
        "Radius": "rG",
        "Albedo": "omega",
        "Barrier": "barr",
        "SurfaceSites": "sites",
        "HOPRatio": "hop",
        "MonoLayers": "nMono",
        "DutyCycle": "duty",
        "CRDesorptionTemperature": "Tcr",
        "BranchRatio": "branch",
    }
    user_var = [
        "double mant = mantles(y)",
        "double gdens = (y[IDX_GRAIN0I]+y[IDX_GRAINM])",
        "double garea = (4*pi*rG*rG) * gdens",
        "double unisites = sites * (4*pi*rG*rG)",
        "double densites = sites * (4*pi*rG*rG) * gdens",
        "double freq = sqrt((2.0*sites*kerg)/((pi*pi)*amu))",
        "double quan = -2.0*(barr/hbar) * sqrt(2.0*amu*kerg)",
        "double layers = mant/(nMono*densites)",
        "double cov = (mant == 0.0) ? 0.0 : fmin(layers/mant, 1.0/mant)",
    ]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)