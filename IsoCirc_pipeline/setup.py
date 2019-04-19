from setuptools import setup, find_packages
import re, io
from isocirc.__init__ import __version__
from isocirc.__init__ import __program__

# __version__ = re.search(
# 	r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]',  # It excludes inline comment too
# 	io.open('isocirc/__init__.py', encoding='utf_8_sig').read()
# ).group(1)

setup(
	name=__program__,
	# packages=['isocirc'],
	packages=find_packages(),
	version='{}'.format(__version__),
	scripts=['bin/bed2exonGtf',
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
			 'bin/minimap2',
			 'bin/trf409.legacylinux64'])
	],
	entry_points={
		'console_scripts': [
		'{}=isocirc.isocirc:main'.format(__program__), 
		'{}_eval=isocirc.eval_with_anno:main'.format(__program__), 
		'{}_comp=isocirc.isocirc_comp:main'.format(__program__),
		'{}_plot=isocirc.isocirc_plot.py:main'.format(__program__),
		]
	},
	install_requires=[
		'biopython',
		'gffutils',
		'mappy',
		'matplotlib',
		'numpy',
		'pandas',
		'pyfaidx',
		'pysam',
	],
	url='https://github.com/yangao07/PARRIS',
	license='MIT',
	author='yan',
	author_email='yangaoucla@gmail.com',
	description='PARRIS: Profiling and Annotating ciRcular RNA with Iso-Seq',
	long_description=open('README.md').read(),
	long_description_content_type="text/markdown",
)
