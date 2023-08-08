from setuptools import setup
import versioneer

requirements = [ "pandas",
                 "netcdf4",
                 "xarray",
                 "pyyaml",
                 "numba"
]

setup(
    name='deltacd',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Channel Depletions Model",
    license="MIT",
    author="California Department of Water Resources",
    author_email='seshadri.rajagopal@water.ca.gov',
    url='https://github.com/cadwrdeltamodeling/deltacd',
    packages=['deltacd'],
    entry_points={
        'console_scripts': [
            'dcd = deltacd.dcd:main',
            'detaw = deltacd.detaw:main',
            'deltacd2dsm2 = deltacd.utils.deltacd2dsm2:main'
        ]
    },
    install_requires=requirements,
    keywords='deltacd',
    classifiers=[
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.7',
    ]
)
