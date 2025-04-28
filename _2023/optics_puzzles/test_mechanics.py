import unittest

import numpy as np

import mechanics

ORIGIN = np.array([0, 0, 0])

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

    def test_prism(self):
        prism1 = mechanics.Prism(2, 2, 1, 1)
        np.testing.assert_almost_equal(prism1.position, ORIGIN)

    def test_cylinder(self):
        cylinder= mechanics.Cylinder(2, 2, 1)
        np.testing.assert_almost_equal(cylinder.position, ORIGIN)

    def test_compond(self):
        body1 = mechanics.Prism(
            1,
            1, 1, 1,
            initial_position=np.array([-0.5, 0, 0])
        )
        body2 = mechanics.Prism(
            1,
            1, 1, 1,
            initial_position=np.array([+0.5, 0, 0])
        )

        prism = mechanics.Compond(body1, body2)
        assert prism.mass == 2
        np.testing.assert_almost_equal(prism.position, ORIGIN)

        total_inertia = mechanics.IPrism(2, 2, 1, 1)
        np.testing.assert_almost_equal(prism.Ibody, total_inertia)

    def test_compond_cylinder(self):
        # shaft and rotor moved around to make a simple gyro, rotation
        # point at origin and center of mass to right.
        # TODO: axis=RIGHT for interia
        shaft = mechanics.Cylinder(
            2,
            2.5, .1,
            initial_position=np.array([1.25, 0, 0])
        )
        rotor = mechanics.Cylinder(
            2,
            .2, 1,
            initial_position=np.array([2.6, 0, 0])
        )

        gyro = mechanics.Compond(shaft, rotor)

        assert gyro.mass == 4
        np.testing.assert_almost_equal(gyro.position, [1.925, 0, 0])
