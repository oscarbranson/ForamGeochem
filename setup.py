from setuptools import setup, find_packages

setup(name='foramgeochem',
      version='0.0.1',
      description='Tools for converting between foraminifera geochemistry and environmental parameters.',
      url='https://github.com/oscarbranson/foramgeochem',
      author='Oscar Branson',
      author_email='oscarbranson@gmail.com',
      license='MIT',
      packages=find_packages(),
      classifiers=['Development Status :: 4 - Beta',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 3',
                   ],
      install_requires=['numpy',
                        # 'pandas',
                        # 'matplotlib',
                        'uncertainties',
                        ],
      package_data={
        'foramgeochem': ['resources/*'],
      },
      zip_safe=True)
