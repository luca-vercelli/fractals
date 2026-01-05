# fractals

Generation of Koch-like fractal curves. 
A formula for calculus of fractal dimension is provided.

## How fractal dimension is calculated

If an object has linear size $s$ and dimension $d$, its *measure* (whatever that means: lenght, area, volume...) is expected to be proportional to $s^d$.

Here, two kinds of curves are supported:

### A. Fractals made up of equal copies of the same figure

If the whole object has linear size $S$, and it is made up of $n$ equal copies of linear size $s$:

$$ S^d = n s^d $$

This equation can be solved algebrically to find $d=\log_{S/s}n$.

Most curves belong or can be reduced to this category: Koch snowflakes, Sierpinski triangle and carpet, Cantor set, ... 

### B. Fractals made up of copies of the same figure, with different sizes

If the whole object has linear size $S$, and it is made up of $n$ different copies whose linear sizes is $s_i$:

$$ S^d = \sum_{i=1}^n s_i^d $$

This equation can be solved *numerically* to find $d$.

We povide an *asymmetrical Koch snowflake* as an example for this category.

## Requirements

The program requires Python 3 and pylab.

If you are using Debian, Ubuntu, and such:

    sudo apt-get install python3-numpy python3-scipy python3-matplotlib

Or, using pip:

    pip install numpy scipy matplotlib

## Run

In *nix environment, type

    python3 ./draw.py

An interface is shown for selection of curves.

