
from os import path

from setuptools import find_packages, setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

# Get package metadata from 'bamliquidator_batch/__about__.py' file
about = {}
with open(path.join(here, 'pipeline_tools', '__about__.py')) as f:
    exec(f.read(), about)

setup(
    name=about['__title__'],

    version=about['__version__'],

    description=about['__summary__'],
    long_description=long_description,

    url=about['__url__'],

    maintainer=about['__maintainer__'],
    maintainer_email=about['__email__'],

    license=about['__license__'],

    packages=find_packages(exclude=['tests', 'tests.*', '*.tests', '*.tests.*']),

    install_requires=['numpy'],

    package_data={
		'pipeline_tools': ['annotation/*'],
	},

    include_package_data=True,

    classifiers=[
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
    ],

    keywords=['bioinformatics', 'pipeline', 'linlab'],
)
