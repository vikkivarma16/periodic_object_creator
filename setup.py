# setup.py

from setuptools import setup, find_packages

setup(
    name='periodic_object_creator',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'sympy',
    ],
    description='A package to creat 3D, 2D and 1D objects with periodic properties.',
    author='vikkivarma16',
    author_email='vikkivarma16@gmail.com',
    url='https://github.com/vikkivarma16/periodic_object_creator',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)

