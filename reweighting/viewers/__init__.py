# -*- coding: utf-8 -*-
# **************************************************************************
# Module to declare viewers
# Find documentation here: https://scipion-em.github.io/docs/docs/developer/creating-a-viewer
# **************************************************************************
from .viewer_weights import ReweightingWeightsViewer
from .viewer_correlation import ReweightingCorrelationViewer

try:
    from xmipp3.viewers import XmippLogLikelihoodViewer
    from reweighting.protocols import ReweightingProtComputeLikelihood
    XmippLogLikelihoodViewer._targets.extend([ReweightingProtComputeLikelihood])
except ImportError:
    pass
