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

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)