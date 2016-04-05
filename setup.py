from setuptools import setup, find_packages


REQUIRES = [
    'numpy==1.11.0',
    'python-dateutil==2.5.2',
    'pytz==2016.3',
    'scipy==0.17.0',
    'six==1.10.0'
]

EXCLUDE_FROM_PACKAGES = []

setup(
    name='cdpvals',
    version='0.0.1',
    description='Combining p-values from dependent tests',
    url='https://github.com/gungorbudak/cdpvals',
    author='Gungor Budak',
    author_email='gbudak@iupui.edu',
    license='MIT',
    install_requires=REQUIRES,
    packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
    include_package_data=False,
    zip_safe=False,
    entry_points={}
)
