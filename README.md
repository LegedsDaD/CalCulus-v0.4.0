# CalCulus-v0.4.0
CalCulus* is a fast, lightweight scientific computing library written in C++ and exposed to Python using pybind11.

README WILL BE UPDATED SOON DOWNLOAD USING
```
pip install https://github.com/LegedsDaD/CalCulus/releases/download/v0.4.0/calculus-0.4.0-cp313-cp313-win_amd64.whl
```
TEST CODE
```python
import math
from calculus import (
    Scalar, Vec3, Solver,
    sin, cos, tan, exp, log, sqrt, pow
)

print("=" * 60)
print("🧮 CalCulus Library – Full Feature Test")
print("=" * 60)

# --------------------------------------------------
# 1. Scalar + Automatic Differentiation
# --------------------------------------------------
print("\n[1] Scalar Auto-Differentiation")

# f(x) = x^3 + 2x
x = Scalar(3.0, 1.0)   # d/dx = 1
f = x * x * x + x * 2

print("f(x) = x^3 + 2x")
print("x =", x.val)
print("f(x) value =", f.val)
print("f'(x) derivative =", f.d)
print("Expected f'(3) =", 3 * 3 * 3 + 2)

# --------------------------------------------------
# 2. Trigonometric Functions
# --------------------------------------------------
print("\n[2] Trigonometric Functions")

x = Scalar(math.pi / 4, 1.0)

s = sin(x)
c = cos(x)

print("sin(x) =", s.val, " derivative =", s.d)
print("cos(x) =", c.val, " derivative =", c.d)

# --------------------------------------------------
# 3. Exponential, Log, Sqrt
# --------------------------------------------------
print("\n[3] exp, log, sqrt")

x = Scalar(2.0, 1.0)

print("exp(x):", exp(x).val, "derivative:", exp(x).d)
print("log(x):", log(x).val, "derivative:", log(x).d)
print("sqrt(x):", sqrt(x).val, "derivative:", sqrt(x).d)

# --------------------------------------------------
# 4. Power Function
# --------------------------------------------------
print("\n[4] Power Function")

x = Scalar(3.0, 1.0)
p = pow(x, 4)

print("x^4 =", p.val)
print("d/dx x^4 =", p.d)
print("Expected derivative =", 4 * (3 ** 3))

# --------------------------------------------------
# 5. Vec3 Operations
# --------------------------------------------------
print("\n[5] Vec3 Operations")

v1 = Vec3(1, 2, 3)
v2 = Vec3(4, 5, 6)

print("v1 · v2 =", v1.dot(v2))
print("|v1| =", v1.magnitude())

# --------------------------------------------------
# 6. Numerical Differentiation
# --------------------------------------------------
print("\n[6] Numerical Differentiation")

def f_num(x):
    return x**3 + 2*x

d = Solver.differentiate(f_num, 3.0)
print("Numerical derivative at x=3:", d)
print("Expected:", 29)

# --------------------------------------------------
# 7. Numerical Integration (Simpson’s Rule)
# --------------------------------------------------
print("\n[7] Numerical Integration")

def f_int(x):
    return x * x

area = Solver.integrate(f_int, 0, 3, 100)
print("∫ x² dx from 0 to 3 =", area)
print("Expected =", 9)

# --------------------------------------------------
# 8. Root Finding (Newton-Raphson)
# --------------------------------------------------
print("\n[8] Root Finding")

def f_root(x):
    return x**2 - 2

root = Solver.root_newton(f_root, 1.0)
print("Root of x² - 2 =", root)
print("Expected ≈", math.sqrt(2))

# --------------------------------------------------
# 9. Mixed Scalar + Python Math
# --------------------------------------------------
print("\n[9] Mixed Expressions")

x = Scalar(1.5, 1.0)
expr = sin(x) * exp(x) + sqrt(x)

print("Expression value =", expr.val)
print("Expression derivative =", expr.d)

print("\n✅ All tests completed successfully.")
print("=" * 60)
```
