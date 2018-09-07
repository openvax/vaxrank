# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from __future__ import (absolute_import,)

import os
import logging
import re

from setuptools import setup

readme_dir = os.path.dirname(__file__)
readme_path = os.path.join(readme_dir, 'README.md')

try:
    with open(readme_path, 'r') as f:
        readme_markdown = f.read()
except:
    logging.warn("Failed to load %s" % readme_path)
    readme_markdown = ""

try:
    import pypandoc
    readme_restructured = pypandoc.convert(readme_markdown, to='rst', format='md')
except:
    readme_restructured = readme_markdown
    logging.warn("Conversion of long_description from MD to RST failed")
    pass


with open('vaxrank/__init__.py', 'r') as f:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        f.read(),
        re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")

if __name__ == '__main__':
    setup(
        name='vaxrank',
        version=version,
        description="Mutant peptide ranking for personalized cancer vaccines",
        author="Alex Rubinsteyn",
        author_email="alex.rubinsteyn@gmail.com",
        url="https://github.com/hammerlab/vaxrank",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            'six',
            'numpy>=1.14.0',
            'pandas',
            'pyensembl>=1.5.0',
            'varcode>=0.5.9',
            'isovar>=0.8.5',
            'mhctools>=1.5.0',
            'roman',
            'jinja2',
            'pdfkit',
            'pypandoc',
            'shellinford>=0.3.4',
            'xlrd',
            'xlsxwriter',
            'xvfbwrapper',
            'future>=0.16.0',  # needed by pylint
            'astropy',
        ],

        long_description=readme_restructured,
        packages=['vaxrank'],
        package_data={'vaxrank': ['templates/*', 'data/*', 'logging.conf']},
        entry_points={
            'console_scripts': [
                'vaxrank = vaxrank.cli:main'
            ]
        }
    )
