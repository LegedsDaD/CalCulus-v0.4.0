from .calculus_core import Scalar, Vec3, Solver, Constants

def _to_scalar(x):
    # If it's already a Scalar, pass it through.
    # If it's a float/int, wrap it. 
    # Note: Wrapping a float creates a Scalar with derivative 0.0.
    # If user wants AD, they must explicitly create Scalar(v, 1.0) in their code.
    return x if isinstance(x, Scalar) else Scalar(x)

# Expose Scalar methods as module-level functions
def sin(x): return _to_scalar(x).sin()
def cos(x): return _to_scalar(x).cos()
def tan(x): return _to_scalar(x).tan()
def exp(x): return _to_scalar(x).exp()
def log(x): return _to_scalar(x).log()
def sqrt(x): return _to_scalar(x).sqrt()
def abs(x): return _to_scalar(x).abs()

def pow(x, n): return _to_scalar(x).pow(n)

__all__ = [
    "Scalar", "Vec3", "Solver", "Constants",
    "sin", "cos", "tan", "exp", "log", "sqrt", "abs", "pow"
]
