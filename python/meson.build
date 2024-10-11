project('startingmeson', 'cpp',
        version: '0.1'
)

py = import('python').find_installation(pure: false)

cpp = meson.get_compiler('cpp')

incdir = include_directories(['../src', '../submodules/cpp-utils/include'])

fftw = dependency('fftw3')

module_path = 'gatdaem1d'

shared_library('gatdaem1d_ext',
        ['../src/gatdaem1d.cpp', '../src/gatdaem1d.h'],
        include_directories: incdir,
        dependencies: fftw,
        install: true)

python_sources = [
    'gatdaem1d/__init__.py',
    ]

py.install_sources(
    python_sources,
    subdir: module_path
)

