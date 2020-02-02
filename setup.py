from os import path

from setuptools import find_packages, setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md")) as f:
    long_description = f.read()

# Get package metadata from 'bamliquidator_batch/__about__.py' file
about = {}
with open(path.join(here, "pipeline_tools", "__about__.py")) as f:
    exec(f.read(), about)

setup(
    name=about["__title__"],
    version=about["__version__"],
    description=about["__summary__"],
    long_description=long_description,
    url=about["__url__"],
    maintainer=about["__maintainer__"],
    maintainer_email=about["__email__"],
    license=about["__license__"],
    packages=find_packages(exclude=["tests", "tests.*", "*.tests", "*.tests.*"]),
    install_requires=["numpy", "scipy"],
    scripts=[
        "pipeline_tools/scripts/enhancerPromoter_waterfall.R",
        "pipeline_tools/scripts/bamPlot_turbo.R",
        "pipeline_tools/scripts/ROSE2_callSuper.R",
        "pipeline_tools/scripts/ROSE2_stitchOpt.R",
        "pipeline_tools/scripts/clusterEnhancer.R",
        "pipeline_tools/scripts/dynamicEnhancer_plot.R",
        "pipeline_tools/scripts/dynamicEnhancer_rank.R",
    ],
    entry_points={
        "console_scripts": [
            "enhancerPromoter=pipeline_tools.tools.enhancerPromoter:main",
            "out_degree_delta=pipeline_tools.tools.out_degree_delta:main",
            "bamPlot_turbo=pipeline_tools.tools.bamPlot_turbo:main",
            "ROSE2_geneMapper=pipeline_tools.tools.ROSE2_geneMapper:main",
            "ROSE2=pipeline_tools.tools.ROSE2_main:main",
            "ROSE2_META=pipeline_tools.tools.ROSE2_META:main",
            "clusterEnhancer=pipeline_tools.tools.clusterEnhancer:main",
            "dynamicEnhancer=pipeline_tools.tools.dynamicEnhancer:main",
        ]
    },
    package_data={"pipeline_tools": ["annotation/*"]},
    include_package_data=True,
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
    keywords=[
        "bioinformatics",
        "pipeline",
        "linlab",
        "enhancerPromoter",
        "out_degree_delta",
        "bamPlot_turbo",
        "ROSE2_geneMapper",
        "ROSE2",
        "ROSE2_META",
        "clusterEnhancer",
        "dynamicEnhancer",
    ],
)
