"""
Methods for converting ra, dec, and redshift data into 3d.
"""

import warnings
import numpy as np
import astropy.units as u
from astropy.units.quantity import Quantity
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM


class ThreeDData:
    """
    Main class for three dimensional data.
    """

    def __init__(
        self,
        longitude: np.ndarray[Quantity],
        latitude: np.ndarray[Quantity],
        redshift: np.ndarray,
        cosmology: FlatLambdaCDM = None,
    ) -> None:
        if not isinstance(longitude[0], Quantity):
            warnings.warn("No unit for longitude value. Assuming degrees.")
            self.longitude = longitude * u.deg
        else:
            self.longitude = longitude

        if not isinstance(latitude[0], Quantity):
            warnings.warn("No unit for longitude value. Assuming degrees.")
            self.latitude = latitude * u.deg
        else:
            self.latitude = latitude
        self.redshift = redshift
        if cosmology is None:
            warnings.warn("No cosmology given. Assuming vanila LCDM. H0=70, Om0=0.3")
            cosmology = FlatLambdaCDM(H0=70, Om0=0.3)
        self.distances = cosmology.luminosity_distance(redshift)

    def to_cartesian(
        self, frame: str = None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Converts the coordinates into 3D cartesian points
        """
        if frame is None:
            frame = "ICRS"
            warnings.warn(
                "No frame was given. Assuming ICRS (equitorial) coordinate system."
            )
        c = SkyCoord(self.longitude, self.latitude, self.distances, frame=frame)
        x = c.cartesian.x.value
        y = c.cartesian.y.value
        z = c.cartesian.z.value
        return x, y, z
