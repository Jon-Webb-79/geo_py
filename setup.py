# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
        name='project_name',
        version='0.1.0',
        description='give short descritpion of package here',
        long_description=readme,
        author='Jonathan A. Webb',
        author_email='webbja123@gmail.com',
        license=license,
        packages=find_packages()
)
