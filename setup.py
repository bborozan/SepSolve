from setuptools import setup

setup_info = dict(
    name="sepsolve",
    version="0.1.1",
    description="Marker gene selection module",
    license="MIT",

    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",

    author="Bartol Borozan",
    author_email="bborozan@mathos.hr",
    url="https://github.com/bborozan/SepSolve",
    packages=["sepsolve"], 
    install_requires=["gurobipy", "pandas", "scipy", "numpy", "scikit-learn"],
)

setup(**setup_info)