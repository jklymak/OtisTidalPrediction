from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='otis_tide_pred',    # This is the name of your PyPI-package.
    version='0.1',                          # Update the version number for new releases
    description='python for http://volkov.oce.orst.edu/tides/otps.html tidal prediction',
    author_email='jklymak@gmail.com',
    author='Jody Klymak',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jklymak/OtisTidalPrediction",
    packages=setuptools.find_packages()
)
