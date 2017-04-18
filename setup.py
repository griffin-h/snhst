from setuptools import setup

setup(name='snhst',
      author=['Curtis McCully'],
      author_email=['cmccully@lco.global'],
      version=0.9,
      packages=['snhst'],
      setup_requires=[],
      install_requires=['numpy', 'astropy', 'matplotlib', 'drizzlepac', 'astroscrappy', 'reproject',
                        'scipy'],
      tests_require=[],
      scripts=[]
      )
