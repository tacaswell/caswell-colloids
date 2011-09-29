# Copyright 2011, Vinothan N. Manoharan, Thomas G. Dimiduk, Rebecca W. Perry,
# Jerome Fung, and Ryan McGorty, Guangnan Meng
#
# This file is part of Holopy.
#
# Holopy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Holopy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Holopy.  If not, see <http://www.gnu.org/licenses/>.

import track
from holopy.process import detrend
from holopy.third_party.tifffile import TIFFfile
import numpy as np


tif = TIFFfile('pacman_bsa_cml2.tif')
i = 0
for page in tif:
    part = track.findparticles(detrend(page.asarray()[0])[..., np.newaxis],
                                    15, .25,.25)
    np.save('particles_{0}'.format(i), part)
    i+=1


outf.close()



