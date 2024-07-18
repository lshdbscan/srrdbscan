import os
import sys


try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst', format='md')
except :
    long_description = open('README.md', encoding='utf-8').read()

try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    sys.stderr.write('Setuptools not found!\n')
    raise


extra_args = ['-std=c++17', '-march=native', '-O3']
extra_link_args = ['-ltbb']

if sys.platform != 'darwin':
    extra_args += ['-fopenmp']
    extra_link_args += ['-fopenmp']
else:
    extra_args += ['-mmacosx-version-min=10.9', '-stdlib=libc++', '-Xclang', '-fopenmp']
    extra_link_args += ['-lomp']
    os.environ['LDFLAGS'] = '-mmacosx-version-min=10.9'

module = Extension(
    'dbscan_srr',
    sources=['python/wrapper/python_wrapper.cpp', 'src/point.cc', 'src/SRR_LSHDBSCAN.cc', 'src/globals.cc', 'src/transformation.cc', 'src/randomGen.cc'], #, 'src/globals.cc', 'src/transformation.cc'],
    library_dirs=["src/"],
    extra_compile_args=extra_args,
    extra_link_args=extra_link_args,
    runtime_library_dirs=['src/'],
    include_dirs=[ 'src/', 'src/include', 'external/pybind11/include', 'libs', 'third_party/HighFive/include', 'third_party/HighFive/include/highfive'])

setup(
    name='dbscan_srr',
    version='0.1',
    author='Anonymous',
    author_email='lshdbscan@gmail.com',
    url='https://github.com/',
    description=
    'DBSCAN )',
    long_description=long_description,
    license='MIT',
    keywords=
    'lsh locality-sensitive hashing clustering dbscan',
    packages=find_packages(),
    include_package_data=True,
    ext_modules=[module])
