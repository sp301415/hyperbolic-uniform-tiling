from mpmath import *
import cairo
from time import sleep, time
from pprint import pprint

mp.dps = 30
eqeps = 10 ** (-mp.dps // 3)
# mp.pretty = True

"""
Parameters
"""
MAX_ITER = 5
P, Q, R = 3, 4, 4
DRAW_SIZE = 2048
STROKE_COLOR = (0, 0, 0)
STROKE_WIDTH = 2

"""
Options
"""
REMOVE_SAME_LINES = True
DRAW_POLYGONS = "ALL"
AUTO_MAX_ITER = False

IMG = cairo.ImageSurface(cairo.FORMAT_ARGB32, DRAW_SIZE, DRAW_SIZE)
CT = cairo.Context(IMG)
STROKE_WIDTH *= 2 / DRAW_SIZE
CT.translate(DRAW_SIZE / 2, DRAW_SIZE / 2)
CT.scale(DRAW_SIZE / 2, -DRAW_SIZE / 2)

CT.set_source_rgba(0, 0, 0, 0)
CT.rectangle(0, 0, 1, 1)
CT.fill()

CT.arc(0, 0, 1, 0, 2 * pi())
CT.stroke_preserve()
CT.set_source_rgb(255, 255, 255)
CT.fill()

CT.set_line_join(cairo.LINE_JOIN_MITER)


def rotate(z, angle):
    return rect(abs(z), arg(z) + angle)


def parg(x):
    if eqeps < arg(x) < pi() - eqeps:
        return arg(x)
    elif -pi() + eqeps < arg(x) < -eqeps:
        return 2 * pi() + arg(x)
    elif arg(x) > pi() - eqeps or arg(x) < -pi() + eqeps:
        return pi()
    else:
        return 0


def ri(z):
    return (float(z.real), float(z.imag))


def draw_z(z):
    CT.set_source_rgb(255, 0, 0)
    CT.arc(*ri(z), 5 * 2 / DRAW_SIZE, 0, 2 * pi())
    CT.stroke_preserve()
    CT.fill()


class HyperbolicLine:
    def __init__(self, x, y):
        """
        Constructs a Hyperbolic Line between two points x and y.
        If x and y forms a line, then __arg is used.
        If x and y forms a circle, then __center and __r is used.
        """
        self.__x = x
        self.__y = y
        self.__is_circular = True

        if almosteq(x, 0, eqeps):
            self.__is_circular = False
            self.__arg = arg(y)
        elif almosteq(y, 0, eqeps):
            self.__is_circular = False
            self.__arg = arg(x)
        elif almosteq((x / y).imag, 0, eqeps):
            self.__is_circular = False
            self.__arg = arg(x)
        else:
            self.__is_circular = True
            z = conj(1 / x)

            w = z - x
            w /= y - x
            c = (x - y) * (w - abs(w) ** 2) / 2j / w.imag - x

            self.__center = -c
            self.__r = abs(c + x)

    def reflect(self, z):
        """
        Refelcts given z with respect to hyperbolic line.
        """
        if not self.__is_circular:
            return rotate(conj(rotate(z, -self.__arg)), self.__arg)

        else:
            return self.__center + self.__r ** 2 * conj(1 / (z - self.__center))

    def end_points(self):
        """
        Returns the endpoint of given line.
        """
        if not self.__is_circular:
            return (exp(mpc(0, self.__arg)), -exp(mpc(0, self.__arg)))
        else:
            d = abs(self.__center)
            a = (self.__r ** 2 - 1 + d ** 2) / (2 * d)
            h = self.__r ** 2 - a ** 2
            pass

    def draw_line(self, color=STROKE_COLOR, stroke_width=STROKE_WIDTH, rotate=0):
        if not self.__is_circular:
            # assert(arg(self.__x) == self.__arg and arg(self.__y) == self.__arg)
            CT.set_line_width(stroke_width)
            CT.set_source_rgb(*color)
            CT.rotate(rotate)
            CT.move_to(*ri(self.__x))
            CT.line_to(*ri(self.__y))
            CT.stroke()
            CT.rotate(-rotate)
        else:
            # return SVG.circle(draw_z_coord(self.__center), DRAW_SIZE/2 * float(self.__r), fill = 'none', stroke = 'blue')
            # assert(almosteq(abs(self.__x - self.__center), self.__r) and almosteq(abs(self.__y - self.__center), self.__r) )
            def orientation(x, y):
                return (
                    det(
                        matrix(
                            [
                                [1, x.real, x.imag],
                                [1, y.real, y.imag],
                                [1, self.__center.real, self.__center.imag],
                            ]
                        )
                    )
                    > 0
                )

            CT.set_line_width(stroke_width)
            CT.set_source_rgb(*color)
            CT.rotate(rotate)
            if orientation(self.__x, self.__y) > 0:
                CT.arc(
                    *ri(self.__center),
                    self.__r,
                    parg(self.__x - self.__center),
                    parg(self.__y - self.__center),
                )
                CT.stroke()
            else:
                CT.arc(
                    *ri(self.__center),
                    self.__r,
                    parg(self.__y - self.__center),
                    parg(self.__x - self.__center),
                )
                CT.stroke()
            CT.rotate(-rotate)


class SchwarzTriangle:
    def __init__(self, p, q, r):
        """
        Defines a triangle with points p, q, r.
        """
        self.p = chop(p)
        self.q = chop(q)
        self.r = chop(r)
        self.pq = HyperbolicLine(self.p, self.q)
        self.qr = HyperbolicLine(self.q, self.r)
        self.rp = HyperbolicLine(self.r, self.p)
        self.line_reflected = {self.pq: False, self.qr: False, self.rp: False}

    def __eq__(self, other):
        if isinstance(other, SchwarzTriangle):
            return (
                almosteq(self.p, other.p, abs_eps=eqeps)
                and almosteq(self.q, other.q, abs_eps=eqeps)
                and almosteq(self.r, other.r, abs_eps=eqeps)
            )

        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def points(self):
        """
        Returns a tuple of points p, q, r.
        """
        return (self.p, self.q, self.r)

    def lines(self):
        """
        Returns a tuple of lines pq, qr, rp.
        """
        return (self.pq, self.qr, self.rp)

    def draw_triangle(self, rotate):
        for l in self.lines():
            if not self.line_reflected[l]:
                l.draw_line(rotate=rotate)

    def undraw_line(self, line):
        if REMOVE_SAME_LINES:
            self.line_reflected[line] = True

    def __str__(self):
        return "{" + ", ".join((str(self.p), str(self.q), str(self.r))) + "}"

    def __repr__(self):
        return str(self)

    def __hash__(self):
        with workdps(5):
            return hash((self.p, self.q, self.r))


"""
Defines the first Schwartz Triangle
"""
p = 0
pq = tanh(
    acosh(
        (cos(pi() / R) + cos(pi() / P) * cos(pi() / Q))
        / (sin(pi() / P) * sin(pi() / Q))
    )
    / 2
)
q = pq
pr = tanh(
    acosh(
        (cos(pi() / Q) + cos(pi() / P) * cos(pi() / R))
        / (sin(pi() / P) * sin(pi() / R))
    )
    / 2
)
r = rect(pr, pi() / P)

base = SchwarzTriangle(p, q, r)

triangles = []

iter = 0
iter_triangles = [base, SchwarzTriangle(base.p, base.rp.reflect(base.q), base.r)]
iter_vertex = ["Q", "R", "Q"]
iter_triangle_is_internal = [False, False]
triangles = iter_triangles[:]

if AUTO_MAX_ITER:
    hyper_pq = 2 * atanh(pq)

    n = 1
    while True:
        if tanh(hyper_pq * n / 2) > 0.99:
            break
        n += 1

    MAX_ITER = n


CALC_TIME = time()
while True:
    if iter == MAX_ITER:
        break
    iter += 1

    assert len(iter_triangles) - sum(iter_triangle_is_internal) + 1 == len(iter_vertex)

    next_triangles = []
    next_vertex = []
    next_triangle_is_internal = []

    offset = 0
    vertex_i = 0
    for i, t in enumerate(iter_triangles):
        next_triangles_temp = [t]
        next_vertex_temp = []
        next_triangle_is_internal_temp = []

        if iter_triangle_is_internal[i]:
            continue

        j = i
        k = i

        front_internal_cnt = 0
        back_internal_cnt = 0

        while True:
            j = (j - 1) % len(iter_triangle_is_internal)

            if iter_triangle_is_internal[j]:
                front_internal_cnt += 1
            else:
                break

        if front_internal_cnt % 2 == 0:
            front_offset = front_internal_cnt // 2
        else:
            front_offset = front_internal_cnt // 2

        while True:
            k = (k + 1) % len(iter_triangle_is_internal)

            if iter_triangle_is_internal[k]:
                back_internal_cnt += 1
            else:
                break

        if back_internal_cnt % 2 == 0:
            back_offset = back_internal_cnt // 2
        else:
            back_offset = back_internal_cnt // 2 + 1

        if iter_vertex[vertex_i] == "P":
            if iter_vertex[vertex_i + 1] == "Q":
                for j in range(P - 1 - front_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.p, tt.q, tt.pq.reflect(tt.r))
                        next_vertex_temp.append("R")
                        nt.undraw_line(nt.pq)
                    else:
                        nt = SchwarzTriangle(tt.p, tt.rp.reflect(tt.q), tt.r)
                        next_vertex_temp.append("Q")
                        nt.undraw_line(nt.rp)
                    next_triangles_temp.append(nt)

                next_triangles_temp.pop(0)
                next_triangles_temp.reverse()
                next_vertex_temp.reverse()

                for j in range(Q - 2 - back_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.qr.reflect(tt.p), tt.q, tt.r)
                        next_vertex_temp.append("P")
                        nt.undraw_line(nt.qr)
                    else:
                        nt = SchwarzTriangle(tt.p, tt.q, tt.pq.reflect(tt.r))
                        next_vertex_temp.append("R")
                        nt.undraw_line(nt.pq)
                    next_triangles_temp.append(nt)

            elif iter_vertex[vertex_i + 1] == "R":
                for j in range(P - 1 - front_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.p, tt.rp.reflect(tt.q), tt.r)
                        next_vertex_temp.append("Q")
                        nt.undraw_line(nt.rp)
                    else:
                        nt = SchwarzTriangle(tt.p, tt.q, tt.pq.reflect(tt.r))
                        next_vertex_temp.append("R")
                        nt.undraw_line(nt.pq)
                    next_triangles_temp.append(nt)

                next_triangles_temp.pop(0)
                next_triangles_temp.reverse()
                next_vertex_temp.reverse()

                for j in range(R - 2 - back_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.qr.reflect(tt.p), tt.q, tt.r)
                        next_vertex_temp.append("P")
                        nt.undraw_line(nt.qr)
                    else:
                        nt = SchwarzTriangle(tt.p, tt.rp.reflect(tt.q), tt.r)
                        next_vertex_temp.append("Q")
                        nt.undraw_line(nt.rp)
                    next_triangles_temp.append(nt)
            else:
                print("AHH THIS SHOULDNT HAPPEN THO")
                raise ValueError

            next_triangle_is_internal_temp = [False] * len(next_triangles_temp)
            next_triangle_is_internal_temp[P - 1 - front_offset - 1] = True

        elif iter_vertex[vertex_i] == "Q":
            if iter_vertex[vertex_i + 1] == "P":
                for j in range(Q - 1 - front_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.p, tt.q, tt.pq.reflect(tt.r))
                        next_vertex_temp.append("R")
                        nt.undraw_line(nt.pq)
                    else:
                        nt = SchwarzTriangle(tt.qr.reflect(tt.p), tt.q, tt.r)
                        next_vertex_temp.append("P")
                        nt.undraw_line(nt.qr)
                    next_triangles_temp.append(nt)
                next_triangles_temp.pop(0)
                next_triangles_temp.reverse()
                next_vertex_temp.reverse()

                for j in range(P - 2 - back_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.p, tt.rp.reflect(tt.q), tt.r)
                        next_vertex_temp.append("Q")
                        nt.undraw_line(nt.rp)
                    else:
                        nt = SchwarzTriangle(tt.p, tt.q, tt.pq.reflect(tt.r))
                        next_vertex_temp.append("R")
                        nt.undraw_line(nt.pq)
                    next_triangles_temp.append(nt)

            elif iter_vertex[vertex_i + 1] == "R":
                for j in range(Q - 1 - front_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.qr.reflect(tt.p), tt.q, tt.r)
                        next_vertex_temp.append("P")
                        nt.undraw_line(nt.qr)
                    else:
                        nt = SchwarzTriangle(tt.p, tt.q, tt.pq.reflect(tt.r))
                        next_vertex_temp.append("R")
                        nt.undraw_line(nt.pq)
                    next_triangles_temp.append(nt)
                next_triangles_temp.pop(0)
                next_triangles_temp.reverse()
                next_vertex_temp.reverse()

                for j in range(R - 2 - back_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.p, tt.rp.reflect(tt.q), tt.r)
                        next_vertex_temp.append("Q")
                        nt.undraw_line(nt.rp)
                    else:
                        nt = SchwarzTriangle(tt.qr.reflect(tt.p), tt.q, tt.r)
                        next_vertex_temp.append("P")
                        nt.undraw_line(nt.qr)
                    next_triangles_temp.append(nt)
            else:
                print("AHH THIS SHOULDNT HAPPEN THO")
                raise ValueError

            next_triangle_is_internal_temp = [False] * len(next_triangles_temp)
            next_triangle_is_internal_temp[Q - 1 - front_offset - 1] = True

        elif iter_vertex[vertex_i] == "R":
            if iter_vertex[vertex_i + 1] == "P":
                for j in range(R - 1 - front_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.p, tt.rp.reflect(tt.q), tt.r)
                        next_vertex_temp.append("Q")
                        nt.undraw_line(nt.rp)
                    else:
                        nt = SchwarzTriangle(tt.qr.reflect(tt.p), tt.q, tt.r)
                        next_vertex_temp.append("P")
                        nt.undraw_line(nt.qr)
                    next_triangles_temp.append(nt)
                next_triangles_temp.pop(0)
                next_triangles_temp.reverse()
                next_vertex_temp.reverse()

                for j in range(P - 2 - back_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.p, tt.q, tt.pq.reflect(tt.r))
                        next_vertex_temp.append("R")
                        nt.undraw_line(nt.pq)
                    else:
                        nt = SchwarzTriangle(tt.p, tt.rp.reflect(tt.q), tt.r)
                        next_vertex_temp.append("Q")
                        nt.undraw_line(nt.rp)
                    next_triangles_temp.append(nt)

            elif iter_vertex[vertex_i + 1] == "Q":
                for j in range(R - 1 - front_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.qr.reflect(tt.p), tt.q, tt.r)
                        next_vertex_temp.append("P")
                        nt.undraw_line(nt.qr)
                    else:
                        nt = SchwarzTriangle(tt.p, tt.rp.reflect(tt.q), tt.r)
                        next_vertex_temp.append("Q")
                        nt.undraw_line(nt.rp)
                    next_triangles_temp.append(nt)
                next_triangles_temp.pop(0)
                next_triangles_temp.reverse()
                next_vertex_temp.reverse()

                for j in range(Q - 2 - back_offset):
                    tt = next_triangles_temp[-1]
                    if j % 2 == 0:
                        nt = SchwarzTriangle(tt.p, tt.q, tt.pq.reflect(tt.r))
                        next_vertex_temp.append("R")
                        nt.undraw_line(nt.pq)
                    else:
                        nt = SchwarzTriangle(tt.qr.reflect(tt.p), tt.q, tt.r)
                        next_vertex_temp.append("P")
                        nt.undraw_line(nt.qr)
                    next_triangles_temp.append(nt)
            else:
                print("AHH THIS SHOULDNT HAPPEN THO")
                raise ValueError

            next_triangle_is_internal_temp = [False] * len(next_triangles_temp)
            next_triangle_is_internal_temp[R - 1 - front_offset - 1] = True

        if vertex_i != 0:
            next_vertex_temp.pop(0)

        next_vertex += next_vertex_temp
        next_triangles += next_triangles_temp
        next_triangle_is_internal += next_triangle_is_internal_temp
        vertex_i += 1

    triangles += next_triangles
    iter_triangles = next_triangles[:]
    iter_vertex = next_vertex[:]
    iter_triangle_is_internal = next_triangle_is_internal[:]

print("MAX_ITER:", MAX_ITER)
print("CALCULATION ELAPSED:", time() - CALC_TIME)

CALC_TIME = time()

print("TOTAL TILES:", len(triangles))
print("OVERLAPPING TILES:", len(triangles) - len(set(triangles)))

if DRAW_POLYGONS == "P":
    for t in triangles:
        for j in range(P):
            t.qr.draw_line(rotate=2 * pi() / P * (j + 1))
elif DRAW_POLYGONS == "Q":
    for t in triangles:
        for j in range(P):
            t.rp.draw_line(rotate=2 * pi() / P * (j + 1))
elif DRAW_POLYGONS == "R":
    for t in triangles:
        for j in range(P):
            t.pq.draw_line(rotate=2 * pi() / P * (j + 1))
elif DRAW_POLYGONS == "ALL":
    for t in triangles:
        for j in range(P):
            t.draw_triangle(rotate=2 * pi() / P * (j + 1))
else:
    print(f'UNDEFINED PARAMETER "{DRAW_POLYGONS}"')
print("DRAWING ELAPSED:", time() - CALC_TIME)

IMG.write_to_png(f"HyperbolicTilingTest.png")
