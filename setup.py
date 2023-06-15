from setuptools import setup, find_packages

setup(
        name='insist-pista',
        version='1.0.27',
        description="Python Image Simulation and Testing Application",
        url='https://github.com/Jack3690/INSIST',
        author='Avinash CK',
        author_email='avinashck90@gmail.com',
        license='BSD 2-clause',
        package_dir={'': 'src'},
        packages=find_packages(where='src'),
        install_requires=['pandas', 'matplotlib', 'astropy', 'photutils',
                          'numpy', 'seaborn', 'opencv-python', 'reproject',
                          ],
        include_package_data=True,
        package_data={'': ['pista/data/*']},
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.9',
        ],
)
