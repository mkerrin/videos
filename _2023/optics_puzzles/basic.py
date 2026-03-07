from manim_imports_ext import *
from _2023.optics_puzzles.objects import *

from scipy.integrate import solve_ivp
from scipy.spatial import transform
from scipy.spatial.transform import Rotation as Rot
import mechanics
import rk4


np.set_printoptions(precision=3, suppress=True)


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
            state = body.step(dt)
            rotation = body.rotation
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


class CylinderInertia(InteractiveScene):

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
        self.frame.reorient(
            -14,  # theta
            70, # phi
            0,  # gamma
        ).move_to([-0.24, 0.12, 0.04])

        shaftH = 3
        shaftR = 1
        shaft = Cylinder(height=shaftH, radius=shaftR, axis=OUT)
        self.add(shaft)

        shaft.move_to([0, 0, shaftH / 2])

        print("CENTER:", shaft.get_center())

        # Create inertia and body that matches cylinder shape
        inertia = mechanics.ICylinder(1, shaftH, shaftR)
        print("INERTIA:", inertia)
        inertia = mechanics.parallel_axis(1, inertia, shaft.get_center())
        print("INERTIA MOVED UP:", inertia)

        # Save points for shaft so that the rotation matrix works correctly
        shaft.save_points()

        # Create rotation after saving points so we don't get confused
        rot = transform.Rotation.from_rotvec(
            np.pi / 4 * np.array([0, 1, 0])
        ).as_matrix()
        # rot = np.identity(3)
        print(rot)

        body = mechanics.Body(
            1, inertia,
            initial_position=ORIGIN,
            initial_rotation=rot,  # np.identity(3)
        )

        # Update the shaft to match body orientation
        def update_shaft(m: Mobject, dt: float):
            m.apply_initial_points_function(
                lambda points: points.dot(body.rotation.T),
                about_point=body.position
            )

        update_shaft(shaft, 0)

        print(body.q)

        i2 = rot.dot(inertia).dot(rot.T)
        print("INERTIA RIGHT:\n", i2)

        print("BODY INERTIA:\n", body.Ibody)


class CylinderInertia2(InteractiveScene):
    # heavy rotor moved to move

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
        self.frame.reorient(
            -14,  # theta
            70, # phi
            0,  # gamma
        ).move_to([-0.24, 0.12, 0.04])

        shaftH = .2
        shaftR = 1
        shaft = Cylinder(height=shaftH, radius=shaftR, axis=OUT)
        # self.add(shaft)

        shaft2 = Cylinder(height=2, radius=0.2, axis=OUT)
        # self.add(shaft2)

        # shaft.move_to([0, 0, shaftH / 2])
        shaft.move_to([0, 0, 1])

        print("CENTER:", shaft.get_center())

        # Create inertia and body that matches cylinder shape
        inertia = mechanics.ICylinder(1, shaftH, shaftR)
        print("INERTIA:", inertia)
        inertia = mechanics.parallel_axis(1, inertia, shaft.get_center())
        print("INERTIA MOVED UP:", inertia)

        shaftG = Group(shaft, shaft2)
        self.add(shaftG)

        # Save points for shaft so that the rotation matrix works correctly
        shaftG.save_points()

        # Create rotation after saving points so we don't get confused
        rot = transform.Rotation.from_rotvec(
            45 * DEGREES * np.array([0, 1, 0])
        ).as_matrix()
        # rot = np.identity(3)
        print(rot)

        omega = np.array([0, 0, 1])

        body = mechanics.Body(
            1, inertia,
            initial_position=ORIGIN,
            initial_rotation=rot,  # np.identity(3)
            initial_omega=rot.dot(omega),
        )
        body.cm = rot.dot(shaft.get_center())

        Larrow = Arrow(ORIGIN, body.L, fill_color=YELLOW)
        self.add(Larrow)
        dotCM = Sphere(radius=0.2, color=RED)
        dotCM.move_to(body.cm)
        self.add(dotCM)

        tail = TracingTail(dotCM, time_traced=60)
        self.add(tail)

        # Update the shaft to match body orientation
        def update_shaft(m: Mobject, dt: float):
            body.step(dt, force=np.array([0, 0, -9.81]))
            cm = body.rotation.dot(body.cm)
            dotCM.move_to(shaft.get_center())
            print(body.q.normalized())
            rot = Rot.from_quat(body.state[18:22], scalar_first=True).as_matrix()
            print(rot)
            m.apply_initial_points_function(
                lambda points: points.dot(body.rotation.T),
                about_point=body.position
            )
        shaftG.add_updater(update_shaft)
        # update_shaft(shaftG, 0)

        print(body.q)

        i2 = rot.dot(inertia).dot(rot.T)
        print("INERTIA RIGHT:\n", i2)

        print("BODY INERTIA:\n", body.Ibody)


class Angular(InteractiveScene):
    # Basic angular momentum vectors and make sure it lines up

    def construct(self):
        axes = ThreeDAxes(
            x_range=(-4, 4),
            y_range=(-2, 2),
            z_range=(-2, 2),
        )
        axes.add_coordinate_labels(font_size=36)
        z_axis_label = axes.get_axis_label(
            "z", axis=axes.get_z_axis(), edge=UP, direction=DR
        )
        self.add(axes, axes.get_axis_labels(), z_axis_label)
        self.frame.reorient(
            -14,  # theta
            70, # phi
            0,  # gamma
        ).move_to([-0.24, 0.12, 0.04])

        h = 0.5
        r = 2
        rotorC = Cylinder(height=h, radius=r, axis=OUT, color=GREY)
        rotorEnd = Disk3D(radius=r, color=RED)
        rotorEnd.move_to([0, 0, h / 2])
        rotorEnd2 = rotorEnd.copy()
        rotorEnd2.move_to([0, 0, - h / 2])

        dot = Sphere(radius=.4)
        dot.move_to(np.array([r, 0, 0]))

        rotor = Group(rotorC, rotorEnd, rotorEnd2, dot)

        self.add(rotor)

        rot = transform.Rotation.from_rotvec(
            90 * DEGREES * np.array([0, 1, 0])
        ).as_matrix()
        # rot = np.identity(3)

        omega = np.array([0, 0, 2])
        omega = omega.dot(rot.T)

        body = mechanics.CylinderFromShape(
            1, rotorC,
            initial_rotation=rot,
            initial_omega=omega
        )

        L = body.getL(body.state)
        O = body.get_omega()
        Larrow = Arrow(
            start=rotor.get_center(), end=L,
            fill_color=YELLOW
        )
        self.add(Larrow)
        Oarrow = Arrow(
            start=rotor.get_center(), end=O,
            fill_color=RED
        )
        self.add(Oarrow)
        print(L, np.linalg.norm(L))
        print(O, np.linalg.norm(O), O * 2)

        rotor.save_points()

        def update_rotor(m: Mobject, dt: float):
            state = body.step(dt, dotcm=dot.get_center())
            r = body.getRotation(state)
            l = body.getL(state)
            o = body.get_omega()
            cm = rotorC.get_center()
            rotor.apply_initial_points_function(
                lambda points: points.dot(r.T),
                about_point=cm
            )
            Larrow.set_points_by_ends(cm, end=l)
            Oarrow.set_points_by_ends(cm, end=o)

        self.wait()
        rotor.add_updater(update_rotor)


class Pendulum(InteractiveScene):

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
        self.frame.reorient(
            -14,  # theta
            70, # phi
            0,  # gamma
        ).move_to([-0.24, 0.12, 0.04])

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

        bodyShaft = mechanics.Cylinder(
            .1, shaft.height, shaft.radius,
            initial_position=shaft.get_center()
        )
        bodyRotor = mechanics.Cylinder(
            .2, rotor.height, rotor.radius,
            initial_position=rotor.get_center()
        )
        print(bodyShaft.position, bodyRotor.position)
        body = mechanics.Compond(
            bodyShaft, bodyRotor,
            # initial_omega=np.array([0., 0., 0.]),
            initial_omega=np.array([10 * np.pi, 0., 0.]),
        )
        print(body.position)

        cmdot = Sphere(radius=.2, color=RED)
        cmdot.move_to(body.cm)
        self.add(cmdot)

        weight = body.mass * np.array([0, 0, -9.81 / 4])

        gyro.save_points()
        def update_gyro(m: Mobject, dt: float):
            if dt == 0:
                return

            state = body.step(dt, force=weight)

            rotation = body.getRotation(body.state)
            m.apply_initial_points_function(
                lambda points: points.dot(rotation.T),
                about_point=ORIGIN,
            )
            cmdot.move_to(body.position.dot(rotation.T))
        gyro.add_updater(update_gyro)

        Larrow = Arrow(start=ORIGIN, end=body.L, fill_color=RED)
        self.add(Larrow)
        Oarrow = Arrow(start=ORIGIN, end=body.omega, fill_color=PURPLE)
        self.add(Oarrow)

        def update_Larrow(m: Mobject, dt: float):
            if dt == 0:
                return
            L = body.L
            Larrow.set_points_by_ends(ORIGIN, L)
        Larrow.add_updater(update_Larrow)
        def update_Oarrow(m: Mobject, dt: float):
            if dt == 0:
                return
            omega = body.omega
            Oarrow.set_points_by_ends(ORIGIN, omega)
        Oarrow.add_updater(update_Oarrow)


class GyroVectors(InteractiveScene):

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
        self.frame.reorient(
            -14,  # theta
            70, # phi
            0,  # gamma
        ).move_to([-0.24, 0.12, 0.04])


class RotationMatrixes(InteractiveScene):
    # Note sure

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
        self.frame.reorient(
            -14,  # theta
            70, # phi
            0,  # gamma
        ).move_to([-0.24, 0.12, 0.04])

        end = np.array([2, 0, 2])
        rot = transform.Rotation.from_rotvec(
            np.pi / 4 * np.array([0, 1, 0])
        )
        print(rot.as_matrix(), rot.as_quat())
        end = end.dot(rot.as_matrix().T)

        arrow = Arrow(ORIGIN, end)
        self.add(arrow)


class GyroAngularVectors(InteractiveScene):
    # Build on Angular and introduce gyro

    def construct(self):
        axes = ThreeDAxes(
            x_range=(-4, 4),
            y_range=(-2, 2),
            z_range=(-2, 2),
        )
        axes.add_coordinate_labels(font_size=36)
        z_axis_label = axes.get_axis_label(
            "z", axis=axes.get_z_axis(), edge=UP, direction=DR
        )
        self.add(axes, axes.get_axis_labels(), z_axis_label)
        self.frame.reorient(
            -14,  # theta
            70, # phi
            0,  # gamma
        ).move_to([-0.24, 0.12, 0.04])

        # Center of mass of gyro
        cm = np.array([2, 0, 0])
        cmdot = Sphere(radius=0.2, color=RED)
        cmdot.move_to(cm)
        self.add(cmdot)

        # Angular momentum of object
        Larrow = Arrow(ORIGIN, [3, 0, 0], fill_color=YELLOW)
        self.add(Larrow)

        # Weight on CM
        weight = np.array([0, 0, -2])
        Warrow = Arrow(cm, cm + weight, fill_color=PURPLE)
        self.add(Warrow)

        # Torque
        torque = np.cross(cm, weight)
        Tarrow = Arrow(ORIGIN, torque, fill_color=RED)
        self.add(Tarrow)

        rot = Rot.from_rotvec(np.pi / 2 * UP)
        inertia = mechanics.ICylinder(1, 0.5, 1)
        rot_matrix = rot.as_matrix()
        interia = rot_matrix.dot(inertia).dot(rot_matrix.T)
        body = mechanics.Body(
            1, inertia,
            # initial_rotation=rot.as_matrix(),
            initial_omega=np.array([3, 0, 0]),
        )
        body.cm = cm

        rotor = Cylinder(height=0.5, radius=1)
        rotor.apply_points_function(
            lambda points: points.dot(body.rotation),
        )
        cmdot.apply_points_function(
            lambda points: points.dot(body.rotation),
            about_point=ORIGIN,
        )
        rotor.move_to(cm)
        self.add(rotor)

        print("L:", body.L)
        print("cm:", body.cm, body.position)
        print("center:", rotor.get_center())

        def update_scene(m: Mobject, dt: float):
            state = body.step(dt, force=weight)
            rotation = body.rotation
            Larrow.set_points_by_ends(ORIGIN, body.L)
            cm = body.cm.dot(rotation)
            cmdot.move_to(cm)
            print(cm, body.L)
        Tarrow.add_updater(update_scene)


class GyroBasic(InteractiveScene):
    # Build on Angular and introduce gyro

    def construct(self):
        axes = ThreeDAxes(
            x_range=(-4, 4),
            y_range=(-2, 2),
            z_range=(-2, 2),
        )
        axes.add_coordinate_labels(font_size=36)
        z_axis_label = axes.get_axis_label(
            "z", axis=axes.get_z_axis(), edge=UP, direction=DR
        )
        self.add(axes, axes.get_axis_labels(), z_axis_label)
        self.frame.reorient(
            -14,  # theta
            70, # phi
            0,  # gamma
        ).move_to([-0.24, 0.12, 0.04])

        h = 0.25
        r = 1
        shaftH = 2 - h / 2
        shaft = Cylinder(height=shaftH, radius=0.1, axis=OUT, color=GREY, opacity=.4)
        shaft.move_to(np.array([0, 0, shaftH / 2]))

        rotorC = Cylinder(height=h, radius=r, axis=OUT, color=GREY)
        rotorC.move_to(np.array([0, 0, 2]))

        rotorEnd = Disk3D(radius=r, color=RED)
        rotorEnd.move_to([0, 0, shaftH + h])
        rotorEnd2 = rotorEnd.copy()
        rotorEnd2.move_to([0, 0, shaftH])

        # rotor = Group(rotorEnd, rotorEnd2, shaft, rotorC)
        rotor = Group(shaft, rotorC)
        print(rotor.get_center(), rotorC.get_center())

        self.add(rotor)

        rot = transform.Rotation.from_rotvec(
            90 * DEGREES * np.array([0, 1, 0])
        ).as_matrix()
        # rot = np.identity(3)

        print(rotor.get_center(), rotorC.get_center())

        omega = np.array([0, 0, 2])
        omega = omega.dot(rot.T)
        print("OMEGA:", omega)

        body = mechanics.CylinderFromShape(
            1, rotorC,
            initial_rotation=rot,
            initial_omega=omega
        )
        body.cm = rotorC.get_center()
        print("CM:", body.cm)
        weight = body.mass * np.array([0, 0, -9.81])
        body.Ibody = mechanics.parallel_axis(body.mass, body.Ibody, body.cm)

        L = body.getL(body.state)
        O = body.get_omega()
        Larrow = Arrow(
            start=rotorC.get_center(), end=rotorC.get_center() + L,
            fill_color=YELLOW
        )
        self.add(Larrow)
        Oarrow = Arrow(
            start=rotorC.get_center(), end=rotorC.get_center() + O,
            fill_color=RED
        )
        self.add(Oarrow)
        print(L, np.linalg.norm(L))
        print(O, 2 * O, np.linalg.norm(O))

        Farrow = Arrow(
            start=rotorC.get_center(), end=rotorC.get_center() + weight,
            fill_color=PURPLE
        )
        self.add(Farrow)

        t = np.cross(rotorC.get_center(), weight)
        Tarrow = Arrow(
            start=ORIGIN, end=ORIGIN + t,
            fill_color=GREEN
        )
        self.add(Tarrow)

        rotor.save_points()
        def update_gyroX(m: Mobject, dt: float):
            print("DT:", dt, body.getL(body.state), rotorC.get_center())
            state = body.step(dt, force=weight)
            r = body.getRotation(state)
            l = body.getL(state)
            o = body.get_omega()
            rotor.apply_initial_points_function(
                lambda points: points.dot(r.T),
                about_point=ORIGIN
            )
            body.cm = rotorC.get_center()
            Farrow.set_points_by_ends(
                rotorC.get_center(), end=rotorC.get_center() + weight
            )
            Larrow.set_points_by_ends(
                ORIGIN, end=l
            )
            Oarrow.set_points_by_ends(
                ORIGIN, end=o
            )
            t = np.cross(rotorC.get_center(), weight)
            Tarrow.set_points_by_ends(
                ORIGIN, end=t
            )
            print("L:", l, "O:", o, "T:", t)
        rotor.add_updater(update_gyroX)
        def update_gyro(m: Mobject, dt: float, Lrot=None):
            if dt == 0:
                return
            print(rotorC.get_center(), Lrot)
            t = np.cross(rotorC.get_center(), weight)
            Lrot += t * dt
            Larrow.set_points_by_ends(
                ORIGIN, end=Lrot
            )
            body.state[15:18] = L
            state = body.step(dt)
            r = body.getRotation(state)
            rotor.apply_initial_points_function(
                lambda points: points.dot(r.T),
                about_point=ORIGIN
            )
            t = np.cross(rotorC.get_center(), weight)
            Tarrow.set_points_by_ends(
                L, end=L + t
            )
        # rotor.add_updater(lambda m, dt: update_gyro(m, dt, L))


class Gyro(InteractiveScene):

    def construct(self):
        axes = ThreeDAxes(
            x_range=(-4, 4),
            y_range=(-2, 2),
            z_range=(-2, 2),
        )
        axes.add_coordinate_labels(font_size=36)
        z_axis_label = axes.get_axis_label(
            "z", axis=axes.get_z_axis(), edge=UP, direction=DR
        )
        self.add(axes, axes.get_axis_labels(), z_axis_label)
        self.frame.reorient(
            -14,  # theta
            70, # phi
            0,  # gamma
        ).move_to([-0.24, 0.12, 0.04])

        shaftH = 2.5
        shaftR = .1
        rotorH = .2
        rotorR = 1

        shaft = Cylinder(height=shaftH, radius=shaftR, axis=OUT, opacity=0.4)
        rotor = Cylinder(height=rotorH, radius=rotorR, axis=OUT, opacity=0.4)

        rotorEnd = Disk3D(radius=rotorR, color=GREY)
        rotorEnd.move_to([0, 0, shaftH / 2])
        rotorEnd2 = rotorEnd.copy()
        rotorEnd2.move_to([0, 0, shaftH / 2 + rotorH])

        rotor.move_to(np.array([0, 0, shaftH / 2 + rotorH / 2]))

        gyro = Group(shaft, rotor, rotorEnd, rotorEnd2)
        gyro.move_to([0, 0, shaftH / 2 + gyro.get_center()[2]])

        # add gravity to rotor
        shaftM = 1
        rotorM = 2

        bodyShaft = mechanics.Cylinder(
            .1, shaft.height, shaft.radius,
            initial_position=shaft.get_center()
        )
        bodyRotor = mechanics.Cylinder(
            .2, rotor.height, rotor.radius,
            initial_position=rotor.get_center()
        )

        rot = transform.Rotation.from_rotvec(
            10 * DEGREES * np.array([0, 1, 0])
        ).as_matrix()
        # rot = np.identity(3)

        omega = np.array([0, 0, 6.15]).dot(rot.T)

        body = mechanics.Compond(
            bodyShaft, bodyRotor,
            initial_rotation=rot,
            initial_omega=omega,  # np.array([0.2, 0., 2.]),
            # initial_omega=np.array([6 * np.pi, 0., 0.]),
        )

        print("OMEGA:", omega)
        print("L:", body.getL(body.state))

        dot = Sphere(radius=0.15, color=RED)
        dot.move_to(body.cm)
        # self.add(dot)

        print("BEFORE:", dot.get_center())

        gyro2 = Group(gyro, dot)
        self.add(gyro2)

        #gyro2.apply_points_function(
        #    lambda points: points.dot(rot.T),
        #    about_point=ORIGIN
        #)
        #dot.apply_points_function(
        #    lambda points: points.dot(rot.T),
        #    about_point=ORIGIN
        #)

        print("NOW:", dot.get_center())

        gyro2.save_points()

        print(body.position)  # center of mass
        weight = body.mass * np.array([0, 0, -9.81 / 4])

        print(dot.get_center(), dot.get_center_of_mass())
        x = body.position.dot(rot.T)
        print(x)
        Farrow = Arrow(
            start=x, end=x + -2 * OUT,
            fill_color=RED, thickness=6
        )
        self.add(Farrow)
        Tarrow = Arrow(
            start=ORIGIN, end=ORIGIN,
            fill_color=PURPLE, thickness=6
        )
        self.add(Tarrow)
        Larrow = Arrow(
            start=ORIGIN, end=body.getL(body.state),
            fill_color=YELLOW, thickness=6
        )
        #Larrow.set_points_by_ends(
        #    dot.get_center(), dot.get_center() + body.getL(body.state)
        #)
        self.add(Larrow)

        def update_gyro(m: Mobject, dt: float):
            if dt == 0:
                return

            state = body.step(dt, force=weight)

            rotation = body.getRotation(body.state)
            m.apply_initial_points_function(
                lambda points: points.dot(rotation.T),
                about_point=ORIGIN,
            )
            print(body.cm, dot.get_center())
            t = np.cross(dot.get_center(), weight)
            Farrow.set_points_by_ends(
                dot.get_center(), dot.get_center() + -2 * OUT
            )
            Tarrow.set_points_by_ends(dot.get_center(), dot.get_center() + t)
            Larrow.set_points_by_ends(
                dot.get_center(), dot.get_center() + body.getL(body.state)
            )
        gyro2.add_updater(update_gyro)


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
        self.frame.reorient(
            -20,  # theta
            70, # phi
            10,  # gamma
        )

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
            state = body.step(dt)
            rotation = body.rotation
            m.apply_initial_points_function(
                lambda points: points.dot(rotation.T),
                about_point=cm,
            )
        handle.add_updater(update_handler)

        Larrow = Arrow(start=cm, end=body.getL(body.state) - cm)
        self.add(Larrow)
        Oarrow = Arrow(start=cm, end=body.get_omega() - cm)
        self.add(Oarrow)

        def update_Larrow(m: Mobject, dt: float):
            if dt == 0:
                return
            L = body.L
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
            rotation = body.rotation
            m.apply_initial_points_function(
                lambda points: points.dot(rotation.T),
                about_point=cm,
            )
        body_axis.add_updater(update_body_axis)


class Rotation2(InteractiveScene):

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
        self.frame.reorient(
            -20,  # theta
            70, # phi
            10,  # gamma
        )

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

        topBody = mechanics.PrismFromShape(
            1, top, initial_position=top.get_center()
        )
        barBody = mechanics.PrismFromShape(
            1, bar, initial_position=bar.get_center()
        )

        # Add physics to this - mass = 1 for each part
        interia = \
            mechanics.IPrism(1, 1, 5, 1, top.get_center()) + \
            mechanics.IPrism(1, 2, 1, 1, bar.get_center())
        cm = (top.get_center() + bar.get_center()) / 2
        print("INTERIA:", interia, cm)

        self.add(handle)
        handle.save_points()

        bodyOld = mechanics.Body(
            2, interia,
            initial_velocity=np.array([0., 0., 0.]),
            # initial_omega=np.array([0, 2, 0]),
            initial_omega=np.array([np.pi, 0.00002, 0]),
        )
        body = mechanics.Compond(
            barBody, topBody,
            initial_velocity=np.array([0., 0., 0.]),
            # initial_omega=np.array([0, 2, 0]),
            initial_omega=np.array([np.pi, 0.00002, 0]),
        )
        print("body.Ibody:", body.Ibody)
        print("boldOld.Ibody:", bodyOld.Ibody)

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
            state = body.step(dt)
            rotation = body.rotation
            m.apply_initial_points_function(
                lambda points: points.dot(rotation.T),
                about_point=cm,
            )
        handle.add_updater(update_handler)

        Larrow = Arrow(start=cm, end=body.getL(body.state) - cm)
        self.add(Larrow)
        Oarrow = Arrow(start=cm, end=body.get_omega() - cm)
        self.add(Oarrow)

        def update_Larrow(m: Mobject, dt: float):
            if dt == 0:
                return
            L = body.L
            Larrow.set_points_by_ends(cm, L + cm)
        Larrow.add_updater(update_Larrow)
        def update_Oarrow(m: Mobject, dt: float):
            if dt == 0:
                return
            omega = body.omega
            Oarrow.set_points_by_ends(cm, omega + cm)
        Oarrow.add_updater(update_Oarrow)

        def update_body_axis(m: Mobject, dt: float):
            if dt == 0:
                return
            rotation = body.rotation
            m.apply_initial_points_function(
                lambda points: points.dot(rotation.T),
                about_point=cm,
            )
        body_axis.add_updater(update_body_axis)
