from distutils.core import setup
setup(
    name = "geomag",
    packages = ["geomag"],
    data_files = [('geomag', ('geomag/WMM.COF',))],
    version = "0.9",
    description = "Magnetic variation/declination",
    author = "Christopher Weiss",
    author_email = "cmweiss@gmail.com",
    url = "http://url.of.geomag/",
    download_url = "http://url.of.geomag/download/geomag-0.9.tgz",
    keywords = ["magnetic", "variation", "declination"],
    classifiers = [
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
        "Topic :: Utilities",
        ],
    long_description = """\
Magnetic variation/declination
------------------------------

Calculates magnetic variation/declination for any latitude/longitude/altitude,
for any date. Uses the NOAA National Geophysical Data Center, epoch 2010 data.
"""
)
