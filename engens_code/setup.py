import setuptools

setuptools.setup(
     name='engens',  
     version='0.0.1',
      scripts=['engens/core/pdbfix.py'],
      data_files=[('', ['engens/core/template_msa.html'])],
     author="kavrakilab&antuneslab",
     author_email="ac121@rice.edu",
     description="Ensemble Genration Package",
     long_description="Ensemble Genration Package",
   long_description_content_type="text/markdown",
     url="https://github.com/AntunesLab/EnGeNs",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GNU General Public License (GPL)",
         "Operating System :: OS Independent",
     ],
 )