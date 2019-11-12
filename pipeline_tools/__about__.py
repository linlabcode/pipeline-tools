"""Central place for package metadata."""

# NOTE: We use __title__ instead of simply __name__ since the latter would
#       interfere with a global variable __name__ denoting object's name.
__title__ = 'pipeline_tools'
__summary__ = 'Pipeline tools'
__url__ = 'https://github.com/linlabcode/pipeline-tools'

# Semantic versioning is used. For more information see:
# https://packaging.python.org/en/latest/distributing/#semantic-versioning-preferred
__version__ = '1.0.0'

__maintainer__ = 'Jost Vrabic Koren'

__email__ = 'jost.koren@bcm.edu'

__license__ = 'MIT'

__all__ = (
    '__title__', '__summary__', '__url__', '__version__', '__maintainer__',  '__license__',
    '__email__',
)
