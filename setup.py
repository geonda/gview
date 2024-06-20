from setuptools import setup

# Package metadata
name = 'gview'
version = '0.0.1'
description = 'tbd'
author = 'A.Geondzhian, N.Davydov'


# Package dependencies
install_requires = [
    'plotly',
    'numpy',
    'siman'
]


# Package setup
setup(
    name=name,
    version=version,
    description=description,
    author=author,
    packages=['gview'],
    install_requires=install_requires,
)
