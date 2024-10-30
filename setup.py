from setuptools import setup

options = {
    "name" : 'SepSolve',
    "version" : '0.1.0',    
    "description" : 'A example Python package',
    "url" : 'https://github.com/bborozan/SepSolve',
    "author" : 'AAAA',
    "author_email" : 'test@aaa.com',
    "license" : 'MIT',
    "packages" : ['sepsolve'],
    "provides" : ['sepsolve'],
    "package_dir": { "sepsolve" : "sepsolve" },


    "install_requires" : [
        "numpy",
        "scipy",
        "gurobipy",
        "pandas",
    ]
}

setup(**options)