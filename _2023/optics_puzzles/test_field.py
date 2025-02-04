import unittest

import numpy as np

from objects import *


class TestFields(unittest.TestCase):

    def test_points_to_particle_info(self):
        particle = ChargedParticle()

        unit_diffs, norms, adjusted_norms = points_to_particle_info(
            particle, np.array([[1, 0, 0], [2, 0, 0]])
        )

        assert np.array_equal(unit_diffs, np.array([[1, 0, 0], [1, 0, 0]]))
        assert np.array_equal(norms, np.array([[1], [2]]))
        assert np.array_equal(adjusted_norms, np.array([[1], [2]]))
