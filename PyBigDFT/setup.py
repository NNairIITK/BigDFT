from distutils.core import setup
import setuptools
#This is a list of files to install, and where
#(relative to the 'root' dir, where setup.py is)
#You could be more specific.
files = ["things/*"]

setup(name = "Py-BigDFT",
    version = "0.1",
    description = "Python module for BigDFT drivers and analysis",
    author = "Luigi Genovese",
    author_email = "luigi.genovese@cea.fr",
    url = "www.bigdft.org",
    license = 'GNU-GPL',
    packages = setuptools.find_packages(),
    install_requires=['Pyfutile','ase','matplotlib','numpy'],
    #'package' package must contain files (see list above)
    #I called the package 'package' thus cleverly confusing the whole issue...
    #This dict maps the package name =to=> directories
    #It says, package *needs* these files.
    #package_data = {'package' : files },
    #'runner' is in the root.
    #scripts = ["runner"],
    long_description = """
Welcome to the PyBigDFT module. This module contains python programs to drive and analyse BigDFT calculations.
The different calculations which are performed with BigDFT code will be also possible thanks to some of the functions 
which are present in this package.
""",
    #
    #This next part it for the Cheese Shop, look a little down the page.
    classifiers = [
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 3 - Alpha',

    # Indicate who your project is intended for
    'Intended Audience :: End Users/Developers',
    'Topic :: Software Development :: Analysis and Running Tools',

    # Pick your license as you wish (should match "license" above)
     'License :: OSI Approved :: GNU General Public License (GPL)',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    ]     
  ) 
