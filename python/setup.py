from setuptools import setup
from distutils.core import Extension

module = Extension('tongrams',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['../include', '../external/emphf', '../external/essentials'],
                    libraries = ['boost_iostreams', 'boost_regex'],
                    extra_compile_args=['-std=c++17'],
                    sources = ['tongrams.cpp'])

setup(ext_modules = [module])
