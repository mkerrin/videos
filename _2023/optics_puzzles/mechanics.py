import numpy as np
import quaternion

from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation
import rk4

ORIGIN = np.array([0., 0., 0.])

def interia_point(mass, r):
    rx = r[0]
    ry = r[1]
    rz = r[2]
    return mass * np.array([
        [+ry * ry + rz * rz, - rx * ry,           - rx * rz],
        [-ry * rx,           + rx * rx + rz * rz, - ry * rz],
        [-rz * rx,           - rz * ry,           + rx * rx + ry * ry],
    ])

def IPrism(mass, w, h, d, point=ORIGIN):
    # x, y, z
    # interia is around the origin
    inertia = np.array([
        [h * h + d * d, 0.0,           0.0],
        [0.0,           w * w + d * d, 0.0],
        [0.0,           0.0,           w * w + h * h]
    ]) * mass / 12
    return parallel_axis(mass, inertia, point)


def ICylinder(mass, h, r, point=ORIGIN):
    # https://en.wikipedia.org/wiki/List_of_moments_of_inertia#List_of_3D_inertia_tensors
    # z is parallel to cylinder
    inertia = mass * np.array([
        [( 3 * r * r + h * h ) / 12, 0.,                       0.],
        [0,                          (3 * r * r + h * h) / 12, 0.],
        [0,                          0.,                       0.5 * r * r],
    ])
    return parallel_axis(mass, inertia, point)


def ICylinder2(mass, h, r, point=ORIGIN):
    # https://en.wikipedia.org/wiki/List_of_moments_of_inertia#List_of_3D_inertia_tensors
    Iparallel = 0.5 * r * r
    Iperp = 1 / 12 * (h * h + r * r)
    inertia = mass * np.array([
        [Iparallel, 0.,    0.],
        [0,         Iperp, 0.],
        [0,         0.,    Iperp],
    ])
    return parallel_axis(mass, inertia, point)


def ISphere(mass, r):
    Ixx = mass * 2 / 5 * r * r
    inertia = mass * np.array([
        [Ixx, 0.,  0.],
        [0.,  Ixx, 0.],
        [0.,  0.,  Ixx],
    ])
    return inertia


def parallel_axis(mass, inertia, point):
    # parallel asis theorem tensor generalization
    norm2 = point.dot(point)
    return inertia + mass * (norm2 * np.identity(3) - np.array([
        [point[0] * point[0], point[0] * point[1], point[0] * point[2]],
        [point[1] * point[0], point[1] * point[1], point[1] * point[2]],
        [point[2] * point[0], point[2] * point[1], point[2] * point[2]],
    ]))


class Body:

    def __init__(
            self, mass, Ibody,
            initial_position=ORIGIN,
            initial_velocity=ORIGIN,
            initial_rotation=np.identity(3),
            initial_omega=ORIGIN,
    ):
        # constants
        self.mass = mass
        self.Ibody = initial_rotation.dot(Ibody).dot(initial_rotation.T)
        self.Ibodyinv = np.linalg.inv(Ibody)

        # state variables
        self.position = initial_position
        self.velocity = initial_velocity  # momentum

        q = Rotation.from_matrix(initial_rotation).as_quat(scalar_first=True)
        self.q = np.quaternion(*q)
        self.L = self.Ibody.dot(initial_omega)

        # calculated quantities
        self.cm = initial_position
        self.rotation = initial_rotation
        self.omega = initial_omega

        # erh
        self.tk = 0
        self.state = self.getState()
        self.dim = self.state.shape

    def getState(self):
        return np.array(
            self.position.tolist() +  # 3
            self.velocity.tolist() +  # 3
            self.rotation.ravel().tolist() +  # 9
            self.L.tolist() +  # 3
            [self.q.w, self.q.x, self.q.y, self.q.z]  # 4
        )

    def getPosition(self, state):
        return state[0:3]

    def getVelocity(self, state):
        return state[3:6]

    def getRotation(self, state):
        return state[6:15].reshape(3, 3)

    def getL(self, state):  # angular momentum
        return state[15:18]

    def getOmega(self, rotation, L):
        Iinv = rotation.dot(self.Ibodyinv).dot(rotation.T)
        return Iinv.dot(L)

    def get_omega(self):
        rotation = self.getRotation(self.state)
        L = self.getL(self.state)
        return self.getOmega(rotation, L)

    def reorient(self, points):
        r = self.rotation
        return points.dot(r.T)

    def step(self, dt: float, torque=ORIGIN, force=ORIGIN, **kwargs):
        if dt == 0:
            return self.getState()

        state = rk4.rk4(
            dydt, self.tk, self.state, dt, body=self,
            torque=torque, force=force, **kwargs
        )
        self.position = state[0:3]
        self.velocity = state[3:6]
        self.rotation = state[6:15].reshape(3, 3)
        self.L = state[15:18]
        self.q = np.quaternion(*state[18:22])

        self.omega = self.getOmega(self.rotation, self.L)

        self.state = state
        self.tk += dt

        return state


class Cylinder(Body):

    def __init__(self,
                 mass,
                 height,
                 radius,
                 axis=np.array([0, 0, 1]),
                 **kwargs):
        inertia = ICylinder(mass, height, radius)
        super().__init__(mass, inertia, **kwargs)


def CylinderFromShape(mass, shape, **kwargs):
    return Cylinder(mass, shape.height, shape.radius, shape.axis, **kwargs)


class Prism(Body):

    def __init__(self,
                 mass,
                 width,
                 height,
                 depth,
                 **kwargs):
        inertia = IPrism(mass, width, height, depth)
        super().__init__(mass, inertia, **kwargs)


def PrismFromShape(mass, shape, **kwargs):
    # TODO get length from sketch method
    bbox = shape.get_bounding_box()
    width = bbox[:, 0][2] - bbox[:, 0][0]
    height = bbox[:, 1][2] - bbox[:, 1][0]
    depth = bbox[:, 2][2] - bbox[:, 2][0]
    return Prism(mass, width, height, depth, **kwargs)


class Sphere(Body):

    def __init__(self, mass, radius, **kwargs):
        inertia = ISphere(mass, radius)
        super().__init__(mass, inertia, **kwargs)


def SphereFromShape(mass, shape, **kwargs):
    return Sphere(mass, shape.radius, **kwargs)


class Compond(Body):

    def __init__(self, *bodies, **kwargs):
        self.bodies = bodies

        mass = sum([b.mass for b in bodies])
        cm = sum([b.mass * b.position for b in bodies]) / mass
        inertia = np.zeros((3, 3))
        for b in bodies:
            inertia += parallel_axis(b.mass, b.Ibody, b.position)

        # About ORIGIN
        # inertia = parallel_axis(mass, inertia, cm)

        super().__init__(
            mass, inertia,
            initial_position=cm,
            **kwargs
        )


def star(a):
    """Turn vector a into 3x3 matrix such that
       >>> np.dot(star(a), v) = np.cross(a, v)
    """
    return np.array([
        [ 0.,   -a[2],  a[1]],
        [ a[2],     0, -a[0]],
        [-a[1],  a[0],  0]
    ])


def dydt(tk, state, body: Body, force=ORIGIN, torque=ORIGIN, **kwargs):
    xdot = body.getVelocity(state)
    xdot2 = np.zeros(3)

    rotation = body.getRotation(state)
    L = body.getL(state)

    q = np.quaternion(*state[18:22])

    # We need to calculate the torque here as the state of rotation
    # changes and we introduce an error in the calculation and it
    # looks like the simulation gains energy
    position = body.cm.dot(rotation.T)
    torque = np.cross(position, force)
    #torque = np.cross(L, np.array([0, 0, -.2]))
    #print(torque)
    # print(body.cm, torque, force)

    omega = body.getOmega(rotation, L)
    rdot = star(omega).dot(rotation)

    oq = np.quaternion(0, *omega)
    qdot = 0.5 * oq * q

    return np.concatenate((
        xdot, xdot2, rdot.ravel(), torque, [q.w, q.x, q.y, q.z]
    ))


if __name__ == "__main__X":
    body = Body(
        1,
        IPrism(1, 1, 3, 1),
        initial_velocity=np.array([1, 0, 0])
    )
    dt = 1.0 / 60
    t_span = np.arange(0, 3, dt)
    state = body.getState()

    sol = state
    for t in t_span:
        sol = rk4.rk4(dydt, t, sol, dt, body=body)
        print(sol)

    #sol = solve_ivp(dydt, (0, 3), state)
    #print(sol.t)
    #print(sol.y)
    #position = sol.y[:, 0:3]
    # print(position)


if __name__ == "__main__":
    body1 = Body(1, IPrism(1, 1, 1, 1), initial_position=np.array([-0.5, 0, 0]))
    #print(body1.position)
    # -2.5 - -1.5
    body2 = Body(1, IPrism(1, 1, 1, 1), initial_position=np.array([0.5, 0, 0]))
    # -1.5 -  1.5
    # print(body2.position)

    print(IPrism(2, 2, 1, 1))

    cm = (body1.position + body2.position) / (body1.mass + body2.mass)
    print(cm)
    # calc inertia around cm
    #Ibody1 = body1.Ibody + body1.mass * np.linalg.norm(body1.position - cm)
    #Ibody2 = body2.Ibody + body2.mass * np.linalg.norm(body2.position - cm)
    b1p = body1.position
    Rmod = body1.position * body1.position
    print("Rmod:", Rmod)
    Rmod = body1.position.dot( body1.position )
    print("Rmod:", Rmod)
    # https://en.wikipedia.org/wiki/Parallel_axis_theorem
    # Tensor generalization
    Ibody1 = body1.Ibody + body1.mass * (Rmod * np.identity(3) - np.array([
        [b1p[0] *  b1p[0], b1p[0] *  b1p[1], b1p[0] *  b1p[2]],
        [b1p[1] *  b1p[0], b1p[1] *  b1p[1], b1p[1] *  b1p[2]],
        [b1p[2] *  b1p[0], b1p[2] *  b1p[1], b1p[2] *  b1p[2]],
    ]))
    print(IPrism(1, 1, 1, 1, body1.position))

    print(Ibody1)
    #print(Ibody2)
    #print(Ibody1 + Ibody2)
