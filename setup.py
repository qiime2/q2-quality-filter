from setuptools import setup, find_packages


# NOTE: If you are testing your plugin with `q2cli` (i.e. the `qiime` command)
# while you are developing it, you'll need to run `qiime dev refresh-cache` to
# see the latest changes to your plugin reflected in the CLI. You'll need to
# run this command anytime you modify your plugin's interface (e.g.
# add/rename/remove a command or its inputs/parameters/outputs).
#
# Another option is to set the environment variable `Q2CLIDEV=1` so that the
# cache is refreshed every time a command is run. This will slow down the CLI
# while developing because refreshing the cache is slow. However, the CLI is
# much faster when a user installs release versions of QIIME 2 and plugins, so
# this slowdown should only be apparent when *developing* a plugin.
#
# This manual refreshing of the `q2cli` cache is necessary because it can't
# detect when changes are made to a plugin's code while under development (the
# plugin's version remains the same across code edits). This manual refreshing
# of the cache should only be necessary while developing a plugin; when users
# install QIIME 2 and your released plugin (i.e. no longer in development), the
# cache will automatically be updated when necessary.
setup(
    name="q2-quality-filter",
    version="2017.1.0.dev0",
    packages=find_packages(),
    # pandas, q2templates and q2-dummy-types are only required for the dummy
    # methods and visualizers provided as examples. Remove these dependencies
    # when you're ready to develop your plugin, and add your own dependencies
    # (if there are any).
    install_requires=['qiime2 == 2017.2.*', 'q2-dummy-types == 2017.2.*',
                      'q2templates == 2017.2.*', 'pandas'],
    author="Daniel McDonald",
    author_email="wasade@gmail.com",
    url="Website for q2-quality-filter",
    # Visit choosealicense.com for some guidance on picking a license.
    # license="BSD-3-Clause",
    description="Basic FASTQ quality filtering",
    entry_points={
        "qiime2.plugins":
        ["q2-quality-filter=q2_quality_filter.plugin_setup:plugin"]
    },
    # If you are creating a visualizer, all template assets must be included in
    # the package source, if you are not using q2templates this can be removed
    package_data={
        "q2_quality_filter": ["assets/index.html"]
    }
)
