import autolens as al
import autolens.plot as aplt
import numpy as np
from matplotlib import pyplot as plt 
from astropy.cosmology import Planck15

Mpc_in_m = 3.08567758 * 10**22  # Mpc in [m]
Day_in_s = 24 * 3600  # day in second
Arcsec2rad = 2 * np.pi / 360 / 3600  # arc second in radian
Light_v = 299792458  # [m/s]

class time_delay(object):
    def __init__(self, lens_galaxies=None, zl=0.5, zs=1.0, cosmos=Planck15):
        self.lens_galaxies = lens_galaxies
        self.zl = zl
        self.zs = zs
        self.cosmos = cosmos

    def deflection_at_position(self,image_x,image_y):
        """
        Quite stupid here, I don't find a better api to call autolens functionlity~
        Use tracer is reasonable!!!!, pix source
        """
        grid_array = np.zeros((1,1,2))
        grid_array[:,:,0] = np.array([[image_y]])
        grid_array[:,:,1] = np.array([[image_x]])
        mygrid = al.Grid2D.manual_native(
            grid=grid_array, pixel_scales=0.05 #pixel_scales can be set arbitrarily here
        )
        defl_map = self.lens_galaxies.deflections_from_grid(mygrid)
        alphay = defl_map.slim_binned[0,0]
        alphax = defl_map.slim_binned[0,1]
        return alphax, alphay

    def ray_shoot(self, image_x, image_y):
        alphax, alphay = self.deflection_at_position(image_x,image_y)
        return image_x - alphax, image_y - alphay
        
    def geometry_term(self, image_x, image_y):
        src_x, src_y = self.ray_shoot(image_x, image_y)
        return  0.5*((image_x - src_x)**2 + (image_y - src_y)**2)

    def potential_term(self, image_x, image_y):
        """
        Quite stupid here, I don't find a better api to call autolens functionlity~
        """
        grid_array = np.zeros((1,1,2))
        grid_array[:,:,0] = np.array([[image_y]])
        grid_array[:,:,1] = np.array([[image_x]])
        mygrid = al.Grid2D.manual_native(
            grid=grid_array, pixel_scales=0.05  #pixel_scales can be set arbitrarily here
        )
        psi_map = self.lens_galaxies.potential_from_grid(mygrid)
        psi = psi_map.slim_binned[0]
        return psi

    def fermat_potential(self, image_x, image_y):
        return self.geometry_term(image_x, image_y) - self.potential_term(image_x, image_y)

    def time_delay_distance(self):
        Dd = self.cosmos.angular_diameter_distance(self.zl).value #in Mpc unit
        Ds = self.cosmos.angular_diameter_distance(self.zs).value
        Dds = self.cosmos.angular_diameter_distance_z1z2(self.zl, self.zs).value
        return (1+self.zl)*Dd*Ds/Dds * Mpc_in_m

    def time_delay_from_fermat_potential(self,image_x, image_y):
        """
        return time-delay in days
        """
        factor = Arcsec2rad**2 / Day_in_s
        return self.time_delay_distance() / Light_v * self.fermat_potential(image_x, image_y) * factor

if __name__ == "__main__":
    sis_lens = al.mp.SphericalIsothermal(
        centre = (0.0, 0.0),
        einstein_radius = 1.5,
    )
    lens_galaxies = al.Galaxy(
        redshift=0.2, mass=sis_lens
    )

    td_obj = time_delay(lens_galaxies=lens_galaxies, zl=0.2, zs=0.7)
    dt = td_obj.time_delay_from_fermat_potential(1.0,1.0)
    assert np.isclose(dt, -34.93059590498458) #check if results is consistent with lenstronomy
    # print(td_obj.geometry_term(1.0,1.0), td_obj.potential_term(1.0,1.0))
    # print(td_obj.fermat_potential(1.0,1.0))
    # print(td_obj.time_delay_from_fermat_potential(1.0,1.0))