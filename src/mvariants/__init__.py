# -*- coding: utf-8 -*-
"""Contains tools for mutation analysis."""

from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound


from .variant_calls import VariantCall
from .caller import VarScan2, Mutect2, OptionHandler, GATK
from .util import parse_vcf2
from .pre_process import (
    SamtoolsmPileupSingleFile,
    GATKPreprocessor,
    GATKLanePostProcessor,
)
from .effect_prediction import VEP
from .signatures import (
    summary__varscan_matched,
    signature_analysis,
    summary__varscan_matched_dump_exons,
)
