import warnings

from .veloxchemlib import AngularMomentumIntegralsDriver, ElectricDipoleIntegralsDriver


def _electric_dipole_integrals_driver_set_origin(
    self, x: float, y: float, z: float
) -> None:
    warnings.warn(
        "Using the function `veloxchem.ElectricDipoleIntegralsDriver.set_origin` instead of the setter property `veloxchem.ElectricDipoleIntegralsDriver.origin` is deprecated and in future versions will stop working.\n",
        category=FutureWarning,
        stacklevel=2,
    )

    self.origin = [x, y, z]


ElectricDipoleIntegralsDriver.set_origin = _electric_dipole_integrals_driver_set_origin


def _angular_momentum_integrals_driver_set_origin(
    self, x: float, y: float, z: float
) -> None:
    warnings.warn(
        "Using the function `veloxchem.AngularMomentumIntegralsDriver.set_origin` instead of the setter property `veloxchem.AngularMomentumIntegralsDriver.origin` is deprecated and in future versions will stop working.\n",
        category=FutureWarning,
        stacklevel=2,
    )

    self.origin = [x, y, z]


AngularMomentumIntegralsDriver.set_origin = (
    _angular_momentum_integrals_driver_set_origin
)
