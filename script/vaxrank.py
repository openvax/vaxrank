#!/usr/bin/env python
#
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

"""
Script to generate vaccine peptide predictions from somatic cancer variants,
patient HLA type, and tumor RNA-seq data.

Example usage:
    vaxrank
        --vcf somatic.vcf \
        --bam rnaseq.bam \
        --min-vaccine-peptide-length 23 \
        --max-vaccine-peptide-length 27 \
        --output-csv vaccine-peptides.csv
"""


from __future__ import print_function, division, absolute_import
import logging

from vaxrank.args import arg_parser
from vaxrank.vaccine_peptides import vaccine_peptides_dataframe_from_args

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    args = arg_parser.parse_args()
    print(args)

    df_vaccine_peptides = vaccine_peptides_dataframe_from_args(args)
    print(df_vaccine_peptides)
    if args.output_csv:
        df_vaccine_peptides.to_csv(args.output_csv)
