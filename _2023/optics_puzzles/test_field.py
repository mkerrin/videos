import unittest

import numpy as np

from objects import *


class TestParticles(unittest.TestCase):

    def test_particle(self):
        particle = ChargedParticle()
        np.testing.assert_array_almost_equal(
            particle.get_bounding_box(),
            np.array([
                [-0.2, -0.2, -0.2],
                [0, 0, 0],
                [0.2, 0.2, 0.2],
            ]),
        )
        # particle.get_bounding_box()[1]
        assert np.array_equal(particle.get_center(), np.array([0., 0., 0.]))
        assert particle.get_charge() == 1.0
        assert particle.get_radius() == 0.2
        assert particle.get_internal_time() == 0.0
        assert particle.time == 0
        np.testing.assert_almost_equal(particle.time_step, 0.03333333)

        particle.scale(2)
        assert particle.get_radius() == 0.4

        assert np.array_equal(
            particle.get_acceleration(), np.array([0, 0, 0])
        )
        delays = [2]
        particle.get_info_from_delays([], delays)
        particle.get_past_acceleration(delays[0])
        particle.get_past_position(delays[0])

    def test_particle_acceleration_zero(self):
        # e_field OscillateOnYOneDField
        # ShowTheEffectsOfOscillatingCharge
        # changes via the oscillation_function(self, time)
        # particle time_step and recent_positions is updated via
        # increment_clock(dt)
        # get_acceleration == [0, 0, 0] but moving
        # velocity = 1m/s
        particle = ChargedParticle()
        particle_group = Group(particle)

        vel = np.array([1., 0., 0.])
        dt = 1 / 30  # time_step
        particle_group.add_updater(lambda m: m.move_to(vel * particle.time))
        particle_group.update(dt)  # calls updaters
        np.testing.assert_almost_equal(
            particle.get_center(),
            np.array([0.033333333, 0., 0.])
        )
        np.testing.assert_almost_equal(particle.time, 0.03333333)
        particle_group.update(dt)
        np.testing.assert_almost_equal(
            particle.get_center(),
            np.array([0.066666666, 0., 0.])
        )

        np.testing.assert_almost_equal(
            particle.get_acceleration(), np.array([0., 0., 0.])
        )
        np.testing.assert_almost_equal(
            particle.get_velocity(), np.array([1., 0, 0])
        )

    def test_particle_acceleration_one(self):
        particle = ChargedParticle()
        particle_group = Group(particle)

        vel = np.array([1., 0., 0.])
        acc = np.array([1., 0., 0.])
        dt = 1 / 30  # time_step
        particle_group.add_updater(lambda m: m.move_to(
            vel * particle.time + 0.5 * acc * particle.time * particle.time
        ))
        particle_group.update(dt)  # calls updaters
        particle_group.update(dt)
        particle_group.update(dt)

        a = particle.get_acceleration()
        np.testing.assert_almost_equal(a, np.array([1, 0, 0]))

    def test_particle_add_force(self):
        # TODO
        pass

    def test_points_to_particle_info(self):
        particle = ChargedParticle()

        unit_diffs, norms, adjusted_norms = points_to_particle_info(
            particle, np.array([
                [1, 0, 0],
                [2, 0, 0],
                [0, 1, 0],
                [0, 2, 0],
                [0, .00001, 0],  # inside radius
            ])
        )

        assert np.array_equal(
            unit_diffs, np.array([
                [1, 0, 0],
                [1, 0, 0],
                [0, 1, 0],
                [0, 1, 0],
                [0, 1, 0],
            ])
        )
        np.testing.assert_array_equal(
            norms, np.array([[1], [2], [1], [2], [0.00001]])
        )
        assert np.array_equal(norms, np.array([
            [1], [2], [1], [2], [0.00001]
        ]))
        np.testing.assert_array_almost_equal(adjusted_norms, np.array([
            [1], [2], [1], [2], [4000.000283]
        ]))
        # assert np.array_equal(adjusted_norms, np.array([
        #     [1], [2], [1], [2], [4000.000283]
        # ]))


class TestFields(unittest.TestCase):

    def test(self):
        assert True
