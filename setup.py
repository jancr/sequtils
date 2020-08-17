import setuptools
import pathlib
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand


with (pathlib.Path(__file__).parent / 'version.txt').open() as f:
    version = f.read().strip()


# add test command to setup.py
class PyTest(TestCommand):
    def run_tests(self):
        import pytest
        errno = pytest.main([])
        sys.exit(errno)
cmdclass = { 'test': PyTest }


# build documentation if sphinx is avalible
command_options = {}
name='sequtils'
try:
    from sphinx.setup_command import BuildDoc
except ImportError:
    pass
else:
    class BuildMan(BuildDoc):
        def initialize_options(self):
            BuildDoc.initialize_options(self)
            self.builder = 'man'

    cmdclass['build_sphinx'] = BuildDoc
    command_options['build_sphinx'] = {'project': ('setup.py', name),
                                       'version': ('setup.py', version)}
    cmdclass['build_man'] = BuildMan
    command_options['build_man'] = {'project': ('setup.py', name),
                                    'version': ('setup.py', version)}


with open("README.rst", "r") as fh:
    long_description = fh.read()
    setuptools.setup(
        name=name,
        version=version,
        scripts=[],
        extras_require={
            'dev': [
                'pytest',
                #  'pytest-pep8',
                'pytest-cov',
                'sphinx-autodoc-typehints',
                'sphinx-rtd-theme',
            ]},
        author="Jan Christian Refsgaard",
        author_email="jancrefsgaard@gmail.com",
        description="Biological Sequence Utility Scripts",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/jancr/sequtils",
        #  packages=setuptools.find_packages(),
        packages=['sequtils'],
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        cmdclass = cmdclass,
        command_options = command_options
    )
