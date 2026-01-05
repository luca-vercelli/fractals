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
from enum import Enum

# all computational weight are is in pylab.plot()

class FractalKind(Enum):
    KOCH = "A"
    NOT_EQ_KOCH = "B"

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
    """
    Position, orientation, dimension of any segment (may represent a 2D object as well)
    This class keeps:
        x, y : cartesian coordinates of initial point
        o_alpha : orientation (angle) of the segment
        o_rho : length of the segment
        alpha, rho : polar coordinates of initial point
    """

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
            self.complexity = None  # no elements
        elif len(sizes) == 1:
            self.complexity = FractalKind.KOCH  # all elements with same size
            self.common_size = sizes.pop()
        else:
            self.complexity = FractalKind.NOT_EQ_KOCH  # elements with different sizes

    def __str__(self):
        return self.description

    def __repr__(self):
        return self.__str__()

    def get_dimension(self):
        """ Calculate fractal dimension """
        if self.complexity == FractalKind.KOCH:
            # here: S**d = n * s**d
            return log(len(self.edges)) / log(self.full_size / self.common_size)
        elif self.complexity == FractalKind.NOT_EQ_KOCH:
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
            return None

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

    def draw(self, num_iter):
        """ Plot this fractal with standard parameters """
        self._draw(0, 0, 0, 1, num_iter)

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
        degree = pi / 180 * 80  # you may try changing this parameter from 0 to 90 degrees
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


class KochAsymm (FractalDrawer):

    def __init__(self):
        FractalDrawer.__init__(self,
            [
                Position(0, 0, 0, 2),
                Position(2, 0, pi/3, 1),
                Position(2.5, sqrt(3)/2, -pi/3, 1),
                Position(3, 0, 0, 1)
            ],
            4,
            "Kock asymmetric snowflake",
            6
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
            "Peano/Hilbert curve (simplified)",
            6
            )


class App(Tk):
    """ call .mainloop() to launch """
    def __init__(self, available_fractals):
        Tk.__init__(self)
        self.available_fractals = available_fractals
        self.selected_fractal = None
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
        self.lb.bind("<<ListboxSelect>>", self.on_select)
        self.lb.place(x=20, y=20)
        self.status = StringVar()
        self.label = Label(self.mainframe, textvariable=self.status)
        self.label.place(x=20, y=230)
        self.input_var = StringVar()
        self.input_entry = Entry(self.mainframe, textvariable=self.input_var)
        self.input_entry.place(x=200, y=100, width=100)
        self.input_button = Button(self.mainframe, text="Draw", command=self.on_click)
        self.input_button.place(x=200, y=200)
        # Ensure all matplotlib windows are closed when main window is closed
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def on_select(self, selected):
        try:
            sel = selected.widget.curselection()
            if not sel:
                return
            self.selected_fractal_index = sel[0]
        except Exception:
            return
        self.selected_fractal = fractal = available_fractals[self.selected_fractal_index]
        # populate entry with suggested iterations for the selected fractal
        try:
            self.input_var.set(str(fractal.suggested_iterations))
        except Exception:
            pass
        self.show_expected()

    def show_expected(self):
        num_iter_str = self.input_var.get()
        try:
            num_iter = int(num_iter_str)
        except Exception:
            self.status.set("Invalid number of iterations: %r" % (num_iter_str,))
            return
        string = ""
        if self.selected_fractal is not None:
            string = "Expected %d calls and %d lines for drawing the figure" % \
                 (self.selected_fractal.get_expected_num_calls(num_iter),
                  self.selected_fractal.get_expected_num_lines(num_iter))
        self.status.set(string)

    def on_click(self):
        num_iter_str = self.input_var.get()
        try:
            num_iter = int(num_iter_str)
        except Exception:
            self.status.set("Invalid number of iterations: %r" % (num_iter_str,))
            return
        fractal = available_fractals[self.selected_fractal_index]
        self.draw(fractal, num_iter)

    def draw(self, fractal, num_iter):
        fig = pylab.figure()
        ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))
        pylab.title("Dim. " + str(fractal.get_dimension()) + " kind " + str(fractal.complexity.value))
        fractal.draw(num_iter)
        try:
            fig.canvas.manager.set_window_title(fractal.description)
        except Exception:
            pass
        plotPoints();
        try:
            pylab.show(block=False)
        except TypeError:
            pylab.show()

    def on_close(self):
        try:
            pylab.close('all')
        except Exception:
            pass
        try:
            self.destroy()
        except Exception:
            pass

##### TEST ############################


if __name__ == "__main__":
    available_fractals = []
    available_fractals.append(KochSnowflake())
    available_fractals.append(KochSnowflake80())
    available_fractals.append(KochAsymm())
    available_fractals.append(SierpinskiTriangle())
    available_fractals.append(SierpinskiCarpet())
    available_fractals.append(CantorSet())
    available_fractals.append(CantorDust())
    available_fractals.append(Dragon())
    available_fractals.append(Hilbert())
    available_fractals = sorted(available_fractals, key=lambda x: x.description)
    App(available_fractals).mainloop()
