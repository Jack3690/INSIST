from setuptools import setup, find_packages

setup(
    name          = 'ipista',
    version       = '1.0.0',    
    description   = 'An astronomical image simulation package',
    url           = 'https://github.com/Jack3690/INSIST',
    author        = 'Avinash CK',
    author_email  = 'avinashck90@gmail.com',
    license       = 'BSD 2-clause', 
    package_dir   = {'':'ipista'},
    packages      = find_packages(where='ipista'),          
    install_requires =['pandas','matplotlib','astropy','photutils',
                      'numpy', 'seaborn','scipy','pyqt5' ],
    include_package_data = True,  
    package_data         = {'': ['ipista/data/*']},
   
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
	  'Programming Language :: Python :: 3.9',
    ],
   
)