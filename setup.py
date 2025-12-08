from setuptools import setup, find_packages

setup(
    name='periodic_object_creator',
    version='0.1.1',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "periodic_object_creator": ["*.so"],   # include all .so files
    },
    install_requires=[
        'sympy',
    ],
    description='A package to create 3D, 2D and 1D objects with periodic properties.',
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

