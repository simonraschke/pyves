import os
import re
import sys
import sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.test import test as TestCommand



class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
        self.extra_link_args = ["static"]



def find_clang():
    out = subprocess.check_output(['clang', '--version'])
    for line in out.decode("utf8").splitlines():
        if "InstalledDir" in line:
            path = os.path.join(line.split(" ")[-1], "clang").replace("\\", "/")
            if platform.system() == "Windows":
                path += ".exe"
            return path

    
def find_clangpp():
    out = subprocess.check_output(['clang++', '--version'])
    for line in out.decode("utf8").splitlines():
        if "InstalledDir" in line:
            path = os.path.join(line.split(" ")[-1], "clang++").replace("\\", "/")
            if platform.system() == "Windows":
                path += ".exe"
            return path



class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(
                    e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(
                    r'version\s*([\d.]+)',
                    out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' +
                      extdir, '-DPYTHON_EXECUTABLE=' + sys.executable]
        print(self.debug)
        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += [
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
            ]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get('CXXFLAGS', ''),
            self.distribution.get_version())

        import pprint
        # print("\nenv")
        # pprint.pprint(env)
        print("\ncmake_args")
        pprint.pprint(cmake_args)
        print("\nbuild_args")
        pprint.pprint(build_args)
        print()

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] +
                              cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] +
                              build_args, cwd=self.build_temp)
        print("\n----------------------------------------------------------------------\nPython tests\n")



class CatchTestCommand(TestCommand):
    """
    A custom test runner to execute both Python unittest tests and C++ Catch-
    lib tests.
    """

    def distutils_dir_name(self, dname):
        """Returns the name of a distutils build directory"""
        dir_name = "{dirname}.{platform}-{version[0]}.{version[1]}"
        return dir_name.format(dirname=dname,
                               platform=sysconfig.get_platform(),
                               version=sys.version_info)

    def run(self):
        # Run Python tests
        super(CatchTestCommand, self).run()
        print("\nPython tests complete, now running C++ tests...\n")
        # Run catch tests
        subprocess.call(
            ['./*_test'],
            cwd=os.path.join(
                'build',
                self.distutils_dir_name('temp')),
            shell=True)



import re
VERSIONFILE="pyves/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))



thelibFolder = os.path.dirname(os.path.realpath(__file__))
requirementPath = thelibFolder + '/requirements.txt'
requires = [] # Examples: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirementPath):
    with open(requirementPath) as f:
        requires = f.read().splitlines()



setup(
    name='pyves',
    version=verstr,
    author='Simon Raschke',
    author_email='',
    url='https://github.com/simonraschke/pyves.git',
    description='vesicle python bindings',
    long_description='',
    packages=find_packages(),
    ext_modules=[CMakeExtension('pyves')],
    cmdclass=dict(build_ext=CMakeBuild, test=CatchTestCommand),
    setup_requires=requires,
    install_requires=requires,
    python_requires='~=3.6',
    zip_safe=False,
)