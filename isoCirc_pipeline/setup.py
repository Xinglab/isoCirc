from setuptools import setup, find_packages
import re, io
from isocirc.__init__ import __version__
from isocirc.__init__ import __program__

setup(
	name=__program__,
	# packages=['isocirc'],
	packages=find_packages(),
        #package_dir = {"": "src"},
	version='{}'.format(__version__),
        scripts=['isocirc/isocirc',
        	'isocirc/isocircPlot',
            'bin/bed2exonGtf',
            'bin/bed2gtf',
            'bin/gtf2bed',
            'bin/gtf2gene',
            'bin/itst_gtf_bed',
            'bin/itst_gtf_gtf',
            'bin/isocirc2bed12'
            ],
        data_files=[
			('bin', 
			 ['bin/bedToGenePred',
			 'bin/fxtools',
			 'bin/genePredToBed',
			 'bin/genePredToGtf',
			 'bin/gtfToGenePred',
			 'bin/lordec-correct',
			 'bin/trf409.legacylinux64'])
	],
	# entry_points={
		# 'console_scripts': [
		# '{}_eval=isocirc.hcBSJ_fullIso:main'.format(__program__), 
		# '{}_comp=isocirc.isocirc_comp:main'.format(__program__),
		# '{}_plot=isocirc.isocirc_plot:main'.format(__program__),
		# ] },
	install_requires=[
		'biopython',
		'gffutils',
		'pyinterval',
		'mappy',
		'matplotlib',
		'numpy',
		'pandas',
		'pyfaidx',
		'pysam',
	],
	url='https://github.com/Xinglab/isoCirc',
	license='GLP',
	author='Yan Gao',
	author_email='yangao07@hit.edu.cn',
	description='isoCirc: computational pipeline to identify high-confidence BSJs and full-length circRNA isoforms from isoCirc long-read data',
	long_description=open('README.md').read(),
	long_description_content_type="text/markdown",
)
