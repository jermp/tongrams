from distutils.core import setup, Extension

module1 = Extension('tongram',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['../', '../emphf'],
                    libraries = ['boost_iostreams', 'boost_system', 'boost_filesystem', 'boost_thread',
                                 'boost_regex', 'boost_chrono', 'boost_date_time', 'boost_atomic', 'pthread'],
                    extra_compile_args=['-std=c++11'],
                    # library_dirs = [''],
                    sources = ['tongrammodule.cpp'])

setup (name = 'tongram',
       version = '1.0',
       description = 'Tongrams Ngrams',
       ext_modules = [module1])
