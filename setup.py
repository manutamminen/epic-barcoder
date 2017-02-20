from distutils.core import setup

setup(
    name='epic-barcoder',
    version='0.1dev',
    packages=['epic-barcoder', ],
    license='MIT',
    install_requires=['pandas'],
    long_description=open('README.md').read(),
)
