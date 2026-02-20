#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <cmath>
#include <string>
#include <vector>

namespace py = pybind11;

/* ================= Scalar (Dual Number for Auto-Diff) ================= */
class Scalar {
public:
    double val; // Value
    double d;   // Derivative (du/dx)

    Scalar(double val = 0.0, double d = 0.0) : val(val), d(d) {}

    // Arithmetic: Operator Overloading with Chain Rule
    // (u + v)' = u' + v'
    Scalar operator+(const Scalar& o) const { return { val + o.val, d + o.d }; }
    Scalar operator+(double o) const { return { val + o, d }; }

    // (u - v)' = u' - v'
    Scalar operator-(const Scalar& o) const { return { val - o.val, d - o.d }; }
    Scalar operator-(double o) const { return { val - o, d }; }

    // (u * v)' = u'v + uv'
    Scalar operator*(const Scalar& o) const {
        return { val * o.val, d * o.val + val * o.d };
    }
    Scalar operator*(double o) const { return { val * o, d * o }; }

    // (u / v)' = (u'v - uv') / v^2
    Scalar operator/(const Scalar& o) const {
        if (std::abs(o.val) < 1e-9) throw std::runtime_error("Division by zero");
        return { val / o.val, (d * o.val - val * o.d) / (o.val * o.val) };
    }

    // Power (u^n)' = n * u^(n-1) * u'
    Scalar pow(double n) const {
        return { std::pow(val, n), n * std::pow(val, n - 1) * d };
    }

    // Negation
    Scalar operator-() const { return { -val, -d }; }

    // Comparison
    bool operator==(const Scalar& o) const { return val == o.val; }
    bool operator!=(const Scalar& o) const { return val != o.val; }
    bool operator<(const Scalar& o) const { return val < o.val; }
    bool operator>(const Scalar& o) const { return val > o.val; }

    // Math Functions (Primitives)
    Scalar sin() const { return { std::sin(val), std::cos(val) * d }; }
    Scalar cos() const { return { std::cos(val), -std::sin(val) * d }; }
    Scalar tan() const {
        double s = std::sin(val);
        double c = std::cos(val);
        if (std::abs(c) < 1e-9) throw std::runtime_error("Tan undefined (cos is 0)");
        return { s / c, (1.0 / (c * c)) * d };
    }

    Scalar exp() const { return { std::exp(val), std::exp(val) * d }; }
    Scalar log() const {
        if (val <= 0) throw std::runtime_error("Log domain error");
        return { std::log(val), (1.0 / val) * d };
    }
    Scalar sqrt() const {
        if (val < 0) throw std::runtime_error("Sqrt of negative number");
        return { std::sqrt(val), (0.5 / std::sqrt(val)) * d };
    }
    Scalar abs() const {
        return { std::abs(val), (val >= 0 ? 1.0 : -1.0) * d };
    }
};

/* ================= Vec3 ================= */
class Vec3 {
public:
    double x, y, z;

    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3& o) const { return { x + o.x, y + o.y, z + o.z }; }
    Vec3 operator-(const Vec3& o) const { return { x - o.x, y - o.y, z - o.z }; }
    Vec3 operator*(double s) const { return { x * s, y * s, z * s }; }
    Vec3 operator/(double s) const { return { x / s, y / s, z / s }; }

    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }

    Vec3 cross(const Vec3& o) const {
        return { y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x };
    }

    double magnitude() const { return std::sqrt(dot(*this)); }

    Vec3 normalize() const {
        double m = magnitude();
        return m == 0 ? Vec3() : (*this) * (1.0 / m);
    }

    double distance(const Vec3& o) const { return (*this - o).magnitude(); }

    // Returns angle in radians between two vectors [0, pi]
    double angle(const Vec3& o) const {
        double m = magnitude() * o.magnitude();
        if (m == 0) return 0.0;
        double d_val = dot(o) / m;
        if (d_val > 1.0) d_val = 1.0;
        if (d_val < -1.0) d_val = -1.0;
        return std::acos(d_val);
    }

    // Reflect vector v around normal n: r = v - 2*(v.n)*n
    Vec3 reflect(const Vec3& n) const {
        return (*this) - n * (2 * dot(n));
    }

    // Project this onto o: proj = (this . o / |o|^2) * o
    Vec3 project(const Vec3& o) const {
        double denom = o.dot(o);
        if (denom == 0) return Vec3();
        return o * (dot(o) / denom);
    }

    std::vector<double> to_list() const {
        return { x, y, z };
    }
};

/* ================= Solver ================= */
class Solver {
public:
    // Simpson's 1/3 Rule Integration
    static double integrate(py::function f, double a, double b, int n) {
        if (n <= 0 || n % 2 != 0) throw std::runtime_error("n must be positive and even");
        double h = (b - a) / n;
        double s = f(a).cast<double>() + f(b).cast<double>();
        for (int i = 1; i < n; ++i) {
            double x = a + h * i;
            s += (i % 2 ? 4 : 2) * f(x).cast<double>();
        }
        return s * h / 3.0;
    }

    // Numerical Differentiation (Central Difference)
    static double differentiate(py::function f, double x, double h_val = 1e-5) {
        return (f(x + h_val).cast<double>() - f(x - h_val).cast<double>()) / (2 * h_val);
    }

    // Newton-Raphson Root Finding
    static double root_newton(py::function f, double x0, double tol = 1e-7, int max_iter = 100) {
        double x = x0;
        for (int i = 0; i < max_iter; ++i) {
            double y = f(x).cast<double>();
            double dy = differentiate(f, x);

            if (std::abs(dy) < 1e-12) throw std::runtime_error("Zero derivative encountered");

            double x_new = x - y / dy;
            if (std::abs(x_new - x) < tol) return x_new;
            x = x_new;
        }
        throw std::runtime_error("Newton-Raphson failed to converge");
    }
};

/* ================= Constants ================= */
namespace Constants {
    constexpr double pi = 3.141592653589793;
    constexpr double e = 2.718281828459045;
    constexpr double phi = 1.618033988749895; // Golden Ratio
    constexpr double gamma = 0.577215664901532; // Euler-Mascheroni

    // Physics
    constexpr double c = 299792458.0;
    constexpr double h = 6.62607015e-34;
    constexpr double Na = 6.02214076e23; // Avogadro
    constexpr double k = 1.380649e-23;
    constexpr double g = 9.80665;
}

/* ================= PYBIND ================= */
PYBIND11_MODULE(calculus_core, m) {
    m.doc() = "CalCulus v0.4.0 C++ Scientific Engine with Auto-Diff";

    py::class_<Scalar>(m, "Scalar")
        .def(py::init<double, double>(), py::arg("val"), py::arg("d") = 0.0)
        .def_readwrite("val", &Scalar::val)
        .def_readwrite("d", &Scalar::d)
        .def("__repr__", [](const Scalar& s) {
        return "<Scalar val=" + std::to_string(s.val) + ", d=" + std::to_string(s.d) + ">";
            })
        .def("__float__", [](const Scalar& s) { return s.val; })
        // Arithmetic
        .def("__add__", [](const Scalar& a, const Scalar& b) { return a + b; })
        .def("__add__", [](const Scalar& a, double b) { return a + b; })
        .def("__radd__", [](const Scalar& a, double b) { return a + b; })
        .def("__sub__", [](const Scalar& a, const Scalar& b) { return a - b; })
        .def("__sub__", [](const Scalar& a, double b) { return a - b; })
        .def("__rsub__", [](const Scalar& a, double b) { return Scalar(b) - a; })
        .def("__mul__", [](const Scalar& a, const Scalar& b) { return a * b; })
        .def("__mul__", [](const Scalar& a, double b) { return a * b; })
        .def("__rmul__", [](const Scalar& a, double b) { return a * b; })
        .def("__truediv__", [](const Scalar& a, const Scalar& b) { return a / b; })
        .def("__truediv__", [](const Scalar& a, double b) { return a / b; })
        .def("__rtruediv__", [](const Scalar& a, double b) { return Scalar(b) / a; })
        .def("__neg__", [](const Scalar& a) { return -a; })
        .def("__pow__", [](const Scalar& a, double b) { return a.pow(b); })
        // Comparison
        .def("__eq__", &Scalar::operator==)
        .def("__lt__", &Scalar::operator<)
        // Math
        .def("sin", &Scalar::sin)
        .def("cos", &Scalar::cos)
        .def("tan", &Scalar::tan)
        .def("exp", &Scalar::exp)
        .def("log", &Scalar::log)
        .def("sqrt", &Scalar::sqrt)
        .def("abs", &Scalar::abs)
        .def("pow", &Scalar::pow);

    py::class_<Vec3>(m, "Vec3")
        .def(py::init<double, double, double>(), py::arg("x") = 0, py::arg("y") = 0, py::arg("z") = 0)
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z)
        .def("__repr__", [](const Vec3& v) {
        return "<Vec3 (" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")>";
            })
        .def("__getitem__", [](const Vec3& v, int i) {
        if (i == 0) return v.x; if (i == 1) return v.y; if (i == 2) return v.z;
        throw py::index_error();
            })
        .def("__setitem__", [](Vec3& v, int i, double d) {
        if (i == 0) v.x = d; else if (i == 1) v.y = d; else if (i == 2) v.z = d;
        else throw py::index_error();
            })
        .def("__add__", &Vec3::operator+)
        .def("__sub__", &Vec3::operator-)
        .def("__mul__", &Vec3::operator*)
        .def("__truediv__", &Vec3::operator/)
        .def("dot", &Vec3::dot)
        .def("cross", &Vec3::cross)
        .def("magnitude", &Vec3::magnitude)
        .def("normalize", &Vec3::normalize)
        .def("distance", &Vec3::distance)
        .def("angle", &Vec3::angle)
        .def("reflect", &Vec3::reflect)
        .def("project", &Vec3::project)
        .def("to_list", &Vec3::to_list);

    py::class_<Solver>(m, "Solver")
        .def_static("integrate", &Solver::integrate, py::arg("f"), py::arg("a"), py::arg("b"), py::arg("n") = 100)
        .def_static("differentiate", &Solver::differentiate, py::arg("f"), py::arg("x"), py::arg("h") = 1e-5)
        .def_static("root_newton", &Solver::root_newton, py::arg("f"), py::arg("x0"), py::arg("tol") = 1e-7, py::arg("max_iter") = 100);

    auto cst = m.def_submodule("Constants");
    cst.attr("pi") = Constants::pi;
    cst.attr("e") = Constants::e;
    cst.attr("phi") = Constants::phi;
    cst.attr("gamma") = Constants::gamma;
    cst.attr("c") = Constants::c;
    cst.attr("h") = Constants::h;
    cst.attr("k") = Constants::k;
    cst.attr("Na") = Constants::Na;
    cst.attr("g") = Constants::g;
}