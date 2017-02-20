from distutils.core import setup

setup(
    name='epic_barcoder',
    version='0.1dev',
    packages=['epic_barcoder', ],
    license='MIT',
    install_requires=['pandas',
                      'epride'],
    long_description=open('README.md').read(),
)
