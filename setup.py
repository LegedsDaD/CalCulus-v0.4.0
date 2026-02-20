from setuptools import setup, Extension
import pybind11
import os

# Define the C++ extension module
ext_modules = [
    Extension(
        "calculus.calculus_core",  # Module name (calculus/calculus_core.so)
        sources=["calculus_core.cpp"],
        include_dirs=[pybind11.get_include()],
        language="c++",
        extra_compile_args=["/O2"] if os.name == "nt" else ["-O3"],
    ),
]

setup(
    name="calculus-cpp",
    version="0.4.0",
    author="LegedsDaD",
    description="High-performance scientific computing with Auto-Differentiation",
    long_description="CalCulus v0.4.0: A C++ powered scientific engine for Python featuring Scalar/Vec3 algebra, Solvers, and forward-mode Automatic Differentiation.",
    ext_modules=ext_modules,
    packages=["calculus"],
    install_requires=["pybind11>=2.6.0"],
    zip_safe=False,
    python_requires=">=3.6",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

