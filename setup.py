# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages


setup(
    name="q2-quality-filter",
    version="2017.2.1",
    packages=find_packages(),
    install_requires=['qiime2 == 2017.2.*', 'q2templates == 2017.2.*',
                      'pandas', 'q2-types == 2017.2.*'],
    author="Daniel McDonald",
    author_email="wasade@gmail.com",
    url="https://github.com/wasade/q2-quality-filter",
    description="Basic FASTQ quality filtering",
    entry_points={
        "qiime2.plugins":
        ["q2-quality-filter=q2_quality_filter.plugin_setup:plugin"]
    },
    package_data={
        "q2_quality_filter": ["assets/index.html"]
    }
)
