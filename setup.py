#
# Copyright (C) 2023 Cai club. All rights reserved.
#
from setuptools import setup

setup(
    name="Hairpin",
    version="0.0.2",
    packages=[
        "hairpin",
    ],
    entry_points={
        "console_scripts": [
            "hmc_extractor = hairpin.hmc_extractor:main",
            "mc_extractor = hairpin.mc_extractor:main",
            "hairpin_cut = hairpin.hairpin_cut:main",
            "run_hairpin = hairpin.run_hairpin:main"
        ]
    },
    long_description=open("README.md").read(),
)
