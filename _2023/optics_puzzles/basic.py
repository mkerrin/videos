from manim_imports_ext import *
from _2023.optics_puzzles.objects import *

from scipy.integrate import solve_ivp
import mechanics
import rk4


class Shapes(InteractiveScene):

    def construct(self):
        axes, plane = get_axes_and_plane()
        self.add(axes, plane)
        self.frame.reorient(
            0,  # theta
            45, # phi
            0,  # gamma
        )

        rod = Cylinder(
            height=2,
            radius=.3,
            axis=RIGHT,
        )
        bar = Cylinder(
            height=1.5,
            radius=.6,
            axis=OUT
        )
        bar.move_to(np.array([-1.2, 0, 0]))

        nut = Group(rod, bar)
        nut.rotate(
            angle=45,
            axis=np.array([1., 0., 0.])
        )
        # self.add(nut)

        objs = [
            Sphere(),
            # Torus(),
            Cylinder(),
            Cone(),
            # Line3D(),
            Disk3D(),
            Square3D(),
            Cube(),
            Prism(),
            VCube(),
            VPrism(),
            Dodecahedron(),
        ]
        for index, obj in enumerate(objs):
            obj.move_to(index * 3 * RIGHT)
        self.add(*objs)


class RotationBasic(InteractiveScene):

    def construct(self):
        axes = ThreeDAxes(
            x_range=(-4, 4),
            y_range=(-4, 4),
            z_range=(-2, 2),
            #width=x_unit * (x_range[1] - x_range[0]),
            #height=y_unit * (y_range[1] - y_range[0]),
            #depth=z_unit * (z_range[1] - z_range[0]),
        )
        axes.add_coordinate_labels(font_size=36)
        z_axis_label = axes.get_axis_label(
            "z", axis=axes.get_z_axis(), edge=UP, direction=DR
        )
        self.add(axes, axes.get_axis_labels(), z_axis_label)

        dot1 = Sphere(radius=.2, color=RED_E)
        dot1.move_to(np.array([-2, 0, 0]))
        dot2 = Sphere(radius=.2, color=RED_A)
        dot2.move_to(np.array([+2, 0, 0]))
        line = Line(dot1.get_center(), dot2.get_center(), color=RED)
        bar1 = Group(dot1, line, dot2)

        inertia1 = mechanics.interia_point(2, dot1.get_center()) \
            + \
            mechanics.interia_point(2, dot2.get_center())

        dot1 = Sphere(radius=.2, color=BLUE_E)
        dot1.move_to(np.array([0, -1, 0]))
        dot2 = Sphere(radius=.2, color=BLUE_A)
        dot2.move_to(np.array([0, +1, 0]))
        line = Line(dot1.get_center(), dot2.get_center(), color=BLUE)
        bar2 = Group(dot1, line, dot2)

        obj = Group(bar1, bar2)
        # obj.rotate(.1, axis=np.array([0, 0, 1]))
        self.add(obj)

        inertia2 = mechanics.interia_point(1, dot1.get_center()) \
            + \
            mechanics.interia_point(2, dot2.get_center())

        print(inertia1)
        print(inertia2)

        obj.save_points()

        body = mechanics.Body(
            6, inertia1 + inertia2, initial_omega=np.array([0, TAU, 0])
        )
        def update_body(m: Mobject, dt: float):
            if dt == 0:
                return
            body.state = rk4.rk4(
                mechanics.dydt, body.tk, body.state, dt, body=body
            )
            body.tk += dt
            rotation = body.getRotation(body.state)
            m.apply_initial_points_function(
                lambda points: points.dot(rotation.T),
                about_point=ORIGIN,
            )
        obj.add_updater(update_body)

        tail1 = TracingTail(
            bar1[2], stroke_width=6, time_traced=30).match_color(bar1[2])
        self.add(tail1)

        tail2 = TracingTail(
            bar2[2], stroke_width=6, time_traced=30).match_color(bar2[2])
        self.add(tail2)


class Gyro(InteractiveScene):

    def construct(self):
        axes = ThreeDAxes(
            x_range=(-4, 4),
            y_range=(-2, 2),
            z_range=(-2, 2),
            #width=x_unit * (x_range[1] - x_range[0]),
            #height=y_unit * (y_range[1] - y_range[0]),
            #depth=z_unit * (z_range[1] - z_range[0]),
        )
        axes.add_coordinate_labels(font_size=36)
        z_axis_label = axes.get_axis_label(
            "z", axis=axes.get_z_axis(), edge=UP, direction=DR
        )
        self.add(axes, axes.get_axis_labels(), z_axis_label)

        shaftH = 2.5
        shaftR = .1
        shaft = Cylinder(height=shaftH, radius=shaftR, axis=RIGHT)
        rotorH = .2
        rotorR = 1
        rotor = Cylinder(height=rotorH, radius=rotorR, axis=RIGHT)
        rotorEnd = Disk3D(radius=rotorR, color=GREY)
        rotorEnd.rotate(np.pi / 2, axis=UP)
        rotorEnd.move_to([shaftH / 2, 0, 0])
        rotorEnd2 = rotorEnd.copy()
        rotorEnd2.move_to([shaftH / 2 + rotorH, 0, 0])
        rotor.move_to(np.array([shaftH / 2 + rotorH / 2, 0, 0]))
        gyro = Group(shaft, rotor, rotorEnd, rotorEnd2)
        gyro.move_to([shaftH / 2 + gyro.get_center()[0], 0, 0])
        self.add(gyro)

        # add gravity to rotor
        shaftM = 1
        rotorM = 2
        cm = (shaftM * shaft.get_center() + rotorM * rotor.get_center()) / \
            (shaftM + rotorM)
        print(cm)
        inertia = mechanics.ICyclinder(shaftM, shaftH, shaftR, cm) + \
            mechanics.ICyclinder(rotorM, rotorH, rotorR, cm)


class Rotation(InteractiveScene):

    def construct(self):
        axes = ThreeDAxes(
            x_range=(-4, 4),
            y_range=(-2, 2),
            z_range=(-2, 2),
            #width=x_unit * (x_range[1] - x_range[0]),
            #height=y_unit * (y_range[1] - y_range[0]),
            #depth=z_unit * (z_range[1] - z_range[0]),
        )
        axes.add_coordinate_labels(font_size=36)
        z_axis_label = axes.get_axis_label(
            "z", axis=axes.get_z_axis(), edge=UP, direction=DR
        )
        self.add(axes, axes.get_axis_labels(), z_axis_label)
        #self.frame.reorient(
        #    0,  # theta
        #    45, # phi
        #    0,  # gamma
        #)

        top = Prism(
            width=1,
            height=5,
            depth=1,
        )
        top.move_to(LEFT)

        bar = Prism(
            width=2,
            height=1,
            depth=1,
        )
        print("BAR:", bar.get_bounding_box(), bar.get_center())
        handle = Group(top, bar)
        print("HANDLE:", handle.get_bounding_box(), handle.get_center())

        # Add physics to this - mass = 1 for each part
        interia = \
            mechanics.IPrism(1, 1, 5, 1, top.get_center()) + \
            mechanics.IPrism(1, 2, 1, 1, bar.get_center())
        cm = (top.get_center() + bar.get_center()) / 2
        print(interia, cm)
        self.add(handle)
        handle.save_points()

        body = mechanics.Body(
            2, interia,
            initial_velocity=np.array([0., 0., 0.]),
            # initial_omega=np.array([0, 2, 0]),
            initial_omega=np.array([np.pi, 0.00002, 0]),
        )

        Larrow = Arrow(start=cm, end=body.getL(body.state) - cm)
        self.add(Larrow)
        Oarrow = Arrow(start=cm, end=body.get_omega() - cm)
        self.add(Oarrow)

        x_line = Line(cm, np.array([4.0, 0, 0]) + cm, color=RED)
        y_line = Line(cm, np.array([0, 4.0, 0]) + cm, color=GREEN)
        z_line = Line(cm, np.array([0, 0, 4.0]) + cm, color=BLUE)
        body_axis = Group(x_line, y_line, z_line)
        self.add(body_axis)
        body_axis.save_points()

        dot = Dot(x_line.get_end(), color=x_line.color)
        self.add(dot)
        dot.add_updater(lambda m: m.move_to(x_line.get_end()))

        tail = TracingTail(dot, time_traced=6)
        self.add(tail)

        # TODO: Note that cm and center do not align, and all movement
        # - translation and rotation need to be done relative to cm

        def update_handler(m: Mobject, dt: float):
            if dt == 0:
                return
            body.state = rk4.rk4(
                mechanics.dydt, body.tk, body.state, dt, body=body
            )
            body.tk += dt
            rotation = body.getRotation(body.state)
            m.apply_initial_points_function(
                lambda points: points.dot(rotation.T),
                about_point=cm,
            )
        handle.add_updater(update_handler)

        def update_Larrow(m: Mobject, dt: float):
            if dt == 0:
                return
            L = body.getL(body.state)
            Larrow.set_points_by_ends(cm, L + cm)
        Larrow.add_updater(update_Larrow)
        def update_Oarrow(m: Mobject, dt: float):
            if dt == 0:
                return
            omega = body.get_omega()
            Oarrow.set_points_by_ends(cm, omega + cm)
        Oarrow.add_updater(update_Oarrow)

        def update_body_axis(m: Mobject, dt: float):
            if dt == 0:
                return
            rotation = body.getRotation(body.state)
            m.apply_initial_points_function(
                lambda points: points.dot(rotation.T),
                about_point=cm,
            )
        body_axis.add_updater(update_body_axis)
