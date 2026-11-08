#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""AutoCDAM: fitting of cationic dummy atom models (CDAM) and metal
Lennard-Jones parameters to target properties.

The drivers ``AutoCDAM`` and ``MetalComplexLJDriver`` are also exported at
the top level of ``veloxchem``; the building blocks live here.
"""

from .cache import CachedCalculator
from .calculator import CalculationPlan, Calculator
from .multisystem import MultiSystemCalculator
from .driver import AutoCDAM
from .evaluator import Evaluator
from .candidate import CDAMCandidate
from .generator import CoordinateGenerator
from .framework import FrameworkSystem, FrameworkBaseline
from .geometry import atom_labels
from .bayesian import BayesianFittingDriver
from .search import ParameterPropertyTable, FittingDriver, SklearnGaussianProcessSurrogate
from .metal_complex import LJScanRecord, MetalComplexLJStudy, QMScanDataset, MetalComplexLJDriver
from .plotting import plot_fit
from .reporting import format_table, samples_table, save_table
from .models import (
    CalculationResult,
    EvaluationResult,
    FitResult,
    FitTarget,
    GeneratedBatch,
    OpenMMSystems,
    ParameterSet,
    ParameterSpec,
    SystemData,
)
