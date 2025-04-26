import unittest

import numpy as np

import mechanics

class TestMechanics(unittest.TestCase):

    def test_prism_inertia(self):
        body1 = mechanics.Body(
            1,
            mechanics.IPrism(1, 1, 1, 1),
            np.array([-0.5, 0, 0])
        )

    def test_parallel_axis_prism_inertia(self):
        body1 = mechanics.Body(
            1,
            mechanics.IPrism(1, 1, 1, 1),
            np.array([-0.5, 0, 0])
        )
        body2 = mechanics.Body(
            1,
            mechanics.IPrism(1, 1, 1, 1),
            np.array([+0.5, 0, 0])
        )

        cm = (body1.mass * body1.position + body2.mass * body2.position) / \
            (body1.mass + body2.mass)
        np.testing.assert_almost_equal(cm, np.array([0, 0, 0]))

        total_inertia = mechanics.IPrism(2, 2, 1, 1)
        total_combined = \
            mechanics.parallel_axis(1, body1.Ibody, body1.position) + \
            mechanics.parallel_axis(1, body2.Ibody, body2.position)
        np.testing.assert_almost_equal(total_combined, total_inertia)
