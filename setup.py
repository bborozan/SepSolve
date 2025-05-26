from setuptools import setup

setup_info = dict(
    name="sepsolve",
    version="0.1.0",
    description="Marker gene selection module",
    license="MIT",

    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",

    author="Bartol Borozan",
    author_email="bborozan@mathos.hr",
    url="https://github.com/bborozan/SepSolve",
    packages=["sepsolve"], 
    install_requires=["wheel", "gurobipy", "pandas", "scipy", "numpy"],
)

setup(**setup_info)