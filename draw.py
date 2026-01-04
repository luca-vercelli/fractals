#!/usr/bin/env python3
# -*- coding: latin1 -*-

from __future__ import division

import pylab
from math import *
from scipy.optimize import fsolve   # fsolve(function,root_estimate,args=(), derivative=None)
try:
    from Tkinter import *
except ImportError:
    from tkinter import *

# all computational weight are is in pylab.plot()


######  FUNCTIONS  ##############################

def to_polar(x, y, cx=0, cy=0):
    """ Convert cartesian coordinates to polar, w.r.t. origin (cx, cy) """
    x = x - cx
    y = y - cy
    rho = abs(hypot(x, y))
    alpha = atan2(y, x)
    return [alpha, rho]


def to_cart(alpha, rho, cx=0, cy=0):
    """ Convert polar coordinates to cartesian, w.r.t. origin (cx, cy) """
    x = rho * cos(alpha)
    y = rho * sin(alpha)
    return [x + cx, y + cy]


def line(x1, y1, x2, y2, *args, **kwargs):
    """ Plot line. VERY slow... """
    # print("line "+str(x1)+" "+str(y1)+" "+str(x2)+" "+str(y2)+" ")
    pylab.plot([x1, x2], [y1, y2], *args, **kwargs)


def linePol(x, y, alpha, rho, *args, **kwargs):
    """ Plot line with polar position. VERY slow... """
    [x2, y2] = to_cart(alpha, rho, x, y)
    line(x, y, x2, y2, *args, **kwargs)

points_x = []
points_y = []


def storePoint(x, y):
    """ Save the point in (points_x,points_y) """
    points_x.append(x)
    points_y.append(y)


def plotPoints(*args, **kwargs):
    """ Plot all the points stored in points_x,points_y """
    pylab.plot(points_x, points_y, ".", markersize=1, *args, **kwargs)


def clearPoints():
    """ Initialize points_x,points_y """
    global points_x, points_y
    points_x = []
    points_y = []


######  CLASSES  ##############################


class Position:
    """ Position and orientation of any 2D-object """

    def __init__(self, x, y, o_alpha, o_rho):
        self.x = x
        self.y = y
        self.o_alpha = o_alpha
        self.o_rho = o_rho
        [self.alpha, self.rho] = to_polar(x, y)

    def __str__(self):
        return "(" + str(self.x) + ";" + str(self.y) + ");" + str(self.o_alpha) + ";" + str(self.o_rho)

    def __repr__(self):
        return self.__str__()


class FractalDrawer (object):
    """
    Abstract class for a generic Koch-style fractal
    We assume that the fractal is built starting from a "basic shape"
    Basic shape is made up of a finite number of "basic lines"
    At each iterations each one of such "basic lines" is replaced with the whole fractal
    There is no limit on position and lenght of basic lines
    However we need a "basis" for the fractal, and we assume this is the segment (0,0)-(full_size,0)
    """

    (EMPTY, KOCH, NOT_EQ_KOCH, MULTI_KOCH) = (0, 1, 2, 3)  # Supported kinds of fractals

    def __init__(self, edges, full_size, description, suggested_iterations):
        """
        Constructor
        @param edges: list<Position> contains the "basic lines" building the fractal
        @param full_size: the "basis" of the fractal
        @param description a short description to be used as window title
        """
        self.description = description
        self.edges = edges
        self.full_size = full_size
        self.suggested_iterations = suggested_iterations
        self.common_size = None
        sizes = set([e.o_rho for e in edges])
        # print "sizes=" + str(sizes) #DEBUG
        if not sizes:
            self.complexity = FractalDrawer.EMPTY  # no elements
        elif len(sizes) == 1:
            self.complexity = FractalDrawer.KOCH  # all elements with same size
            self.common_size = sizes.pop()
        else:
            self.complexity = FractalDrawer.NOT_EQ_KOCH  # elements with different sizes

    def __str__(self):
        return self.description

    def __repr__(self):
        return self.__str__()

    def get_dimension(self):
        """ Calculate fractal dimension """
        if self.complexity == FractalDrawer.KOCH:
            # here: S**d = n * s**d
            return log(len(self.edges)) / log(self.full_size / self.common_size)
        elif self.complexity == FractalDrawer.NOT_EQ_KOCH:
            # here: S**d = Sum_i ( s_i**d )
            def func(d):
                ret = self.full_size ** d
                for e in self.edges:
                    ret -= e.o_rho ** d
                return ret
            res = fsolve(func, 1.0)
            try:
                return float(res[0])
            except Exception:
                return float(res)
        else:
            return -1

    def _draw(self, x, y, o_alpha, o_rho, itnum):
        """
        Plot this fractal
        @param x,y,o_alpha,o_rho   position, rotation,dimension of given object
        @param itnum maximum degree of recursion
        """
        # print("qui "+str(x)+" "+str(y)+" "+str(o_alpha)+" "+str(o_rho)+" ")
        if itnum == 0:
            linePol(x, y, o_alpha, o_rho, color='black')
            # storePoint(x, y)
        else:
            scale = o_rho / self.full_size
            # print "qua " + str(scale)
            for e in self.edges:
                o_alpha1 = o_alpha + e.o_alpha
                o_rho1 = e.o_rho * scale
                [x1, y1] = to_cart(o_alpha + e.alpha, e.rho * scale, x, y)
                self._draw(x1, y1, o_alpha1, o_rho1, itnum - 1)

    def draw(self):
        """ Plot this fractal with standard parameters """
        self._draw(0, 0, 0, 1, self.suggested_iterations)

    def get_expected_num_calls(self, itnum=None):
        if itnum is None:
            itnum = self.suggested_iterations
        n = len(self.edges)
        return sum([n ** i for i in range(0, itnum + 1)])

    def get_expected_num_lines(self, itnum=None):
        if itnum is None:
            itnum = self.suggested_iterations
        return len(self.edges) ** itnum


class KochSnowflake (FractalDrawer):

    def __init__(self):
        FractalDrawer.__init__(self,
            [
                Position(0, 0, 0, 1),
                Position(1, 0, pi/3, 1),
                Position(1.5, sqrt(3/4), -pi/3, 1),
                Position(2, 0, 0, 1)
            ],
            3,
            "Kock snowflake",
            5
            )


class KochSnowflake80 (FractalDrawer) :

    def __init__(self):
        degree = pi/180*80  #try changing this parameter from 0 to 90 degrees
        sin_degree = sin(degree)
        cos_degree = cos(degree)
        FractalDrawer.__init__(self,
            [
                Position(0, 0, 0, 1),
                Position(1, 0, degree, 1),
                Position(1 + cos_degree, sin_degree, -degree, 1),
                Position(1 + 2 * cos_degree, 0, 0, 1)
            ],
            2+2*cos_degree,
            "Kock snowflake 80 deg.",
            5
            )


class SierpinskiTriangle (FractalDrawer):

    def __init__(self):
        FractalDrawer.__init__(self,
            [
                Position(0.0, 0.0, 0.0, 1.0),
                Position(0.5, sqrt(3/4), 0.0, 1.0),
                Position(1.0, 0.0, 0.0, 1.0)
            ],
            2,
            "Sierpinski triangle",
            7
            )


class SierpinskiCarpet(FractalDrawer):

    def __init__(self):
        FractalDrawer.__init__(self,
            [
                Position(0.0, 0.0, 0.0, 1.0),
                Position(1.0, 0.0, 0.0, 1.0),
                Position(2.0, 0.0, 0.0, 1.0),
                Position(0.0, 1.0, 0.0, 1.0),
                Position(2.0, 1.0, 0.0, 1.0),
                Position(0.0, 2.0, 0.0, 1.0),
                Position(1.0, 2.0, 0.0, 1.0),
                Position(2.0, 2.0, 0.0, 1.0),
            ],
            3,
            "Sierpinski carpet (slow)",
            5
            )


class CantorSet (FractalDrawer):
    def __init__(self):
        FractalDrawer.__init__(self,
            [
                Position(0, 0, 0, 1),
                Position(2, 0, 0, 1)
            ],
            3,
            "Cantor set",
            6
            )


class CantorDust(FractalDrawer):
    def __init__(self):
        FractalDrawer.__init__(self,
            [
                Position(0, 0, 0, 1),
                Position(2, 0, 0, 1),
                Position(0, 2, 0, 1),
                Position(2, 2, 0, 1)
            ],
            3,
            "Cantor dust",
            6
            )


class Dragon(FractalDrawer):
    def __init__(self):
        FractalDrawer.__init__(self,
            [
                Position(0, 0, -pi/4, sqrt(2)),
                Position(1, -1, pi/4, 2*sqrt(2)),
                Position(3, 1, -pi/4, sqrt(2))
            ],
            4,
            "Dragon",
            7
            )


class Peano(FractalDrawer):
    def __init__(self):
        FractalDrawer.__init__(self,
            [
                Position(0, 0, 0, 1),
                Position(1, 0, pi/2, 1),
                Position(1, 1, pi, 1),
                Position(0, 1, pi/2, 1),
                Position(0, 2, 0, 3),
                Position(3, 2, -pi/2, 1),
                Position(3, 1, pi, 1),
                Position(2, 1, -pi/2, 1),
                Position(2, 0, 0, 1)
                # don't work
            ],
            3,
            "Peano curve",
            2
            )


class Hilbert(FractalDrawer):
    def __init__(self):
        FractalDrawer.__init__(self,
            [
                Position(0, 1, -pi/2, 1),
                Position(0, 1, 0, 1),
                Position(1, 1, 0, 1),
                Position(2, 0, pi/2, 1)
#                Position(2, 0, -pi/2, 1),
#                Position(1, -1, 0, 1),
#                Position(0, -1, 0, 1),
#                Position(0, -1, pi/2, 1)
            ],
            2,
            "Hilbert curve (simplified)",
            6
            )


class App(Tk):
    """ call .mainloop() to launch """
    def __init__(self, available_fractals):
        Tk.__init__(self)
        self.available_fractals = available_fractals
        self.selected_fractal_index = None
        self.initUI()

    def initUI(self):
        self.geometry("400x250+300+300")
        self.title("Choose fractal")
        self.mainframe = Frame(self, background="gray")
        self.mainframe.pack(fill=BOTH, expand=1)
        self.lb = Listbox(self.mainframe)
        for x in self.available_fractals:
            self.lb.insert(END, x)
        self.lb.bind("<<ListboxSelect>>", self.onSelect)
        self.lb.place(x=20, y=20)
        self.status = StringVar()
        self.label = Label(self.mainframe, textvariable=self.status)
        self.label.place(x=20, y=220)

    def onSelect(self, selected):
        self.selected_fractal_index = selected.widget.curselection()[0]
        fractal = available_fractals[self.selected_fractal_index]
        self.draw(fractal)

    def draw(self, fractal):
        fig = pylab.figure()
        ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
        string = "Expected %d calls and %d lines for drawing the figure" % \
                 (fractal.get_expected_num_calls(), fractal.get_expected_num_lines())
        print(string)
        self.status.set(string)  # FIXME need refresh...
        pylab.title(fractal.get_dimension())
        fractal.draw()
        try:
            fig.canvas.manager.set_window_title(fractal.description)
        except Exception:
            pass
        plotPoints();
        pylab.show()

##### TEST ############################


if __name__ == "__main__":
    available_fractals = []
    available_fractals.append(KochSnowflake())
    available_fractals.append(KochSnowflake80())
    available_fractals.append(SierpinskiTriangle())
    available_fractals.append(SierpinskiCarpet())
    available_fractals.append(CantorSet())
    available_fractals.append(CantorDust())
    available_fractals.append(Dragon())
    #available_fractals.append(Peano())
    available_fractals.append(Hilbert())
    #TODO: Fern; Peano
    available_fractals = sorted(available_fractals, key=lambda x: x.description)
    App(available_fractals).mainloop()
