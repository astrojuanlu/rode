# `rode`

Experiments on doing various mathematical calculations with Python and Rust.

## Euler method

```
$ uv run python euler.py
...
<class 'numpy.ndarray'>
[1.         0.93653348 0.87709496 0.8214288  0.76929557 0.72047106
 0.67474527 0.63192154 0.59181568 ...
```

## Farnocchia orbit propagation algorithm

[TBC]

## Plate deflection

$$
w(x, y; \xi, \eta) = \frac{4 P_c}{\pi^4 D L_x L_y} \sum_{m=1}^\infty \sum_{n=1}^\infty \frac{\sin{\frac{m \pi x}{L_x}} \sin \frac{n \pi y}{L_y} \sin{\frac{m \pi \xi}{L_x}} \sin \frac{n \pi \eta}{L_y}}{\left(\left(\frac{m}{L_x}\right)^2 + \left(\frac{n}{L_y}\right)^2\right)^2}
$$

```
$ uv run python plate.py
...
```
