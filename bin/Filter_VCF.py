#!/usr/bin/env python
"""Apply filters to a VCF on a per-site basis. Edit the constants defined near
the top of the script to adjust the filtering thresholds. Will also remove any
length polymorphism (indels). Takes one argument:
    1) VCF to filter (gzipped)"""

import sys
import gzip

try:
    vcf_in = sys.argv[1]
except IndexError:
    sys.stderr.write(__doc__ + '\n')
    sys.exit(1)

# The filering constants
MIN_VAR_QUAL = 30
MIN_GENO_DP = 10
MAX_SITE_MISSING = 0.1


def geno_depth_flt(vcf_line, mindp):
    """Return a modified VCF record with the low-depth genotypes replaced with
    missing calls."""
    # Keep the first 9 elements of the VCF record the same
    new_vcf = vcf_line[:9]
    # Calcuate the genotype columns and the depth columns
    fmt_spec = vcf_line[8].split(':')
    gt_col = fmt_spec.index('GT')
    dp_col = fmt_spec.index('DP')
    # Iterate through the per-sample information, checking the genotype and
    # depth information.
    for samp in vcf_line[9:]:
        # Freebayes codes missing data as a single dot (.), as opposed to a
        # colon-separated list of dots. If we see a missing call, we will just
        # leave it as-is.
        if samp == '.':
            new_vcf.append(samp)
        else:
            s_dat = samp.split(':')
            try:
                dp = int(s_dat[dp_col])
            except ValueError:
                # If for some reason the depth is not a number, we will just
                # set it to 0.
                dp = 0
            if dp < mindp:
                # We will use a single dot (.) as a missing genotype call.
                new_gt = '.'
            else:
                # If the site has sufficient depth, we keep the genotype call
                # as-is.
                new_gt = s_dat[gt_col]
            # And overwrite the genotype call with the filtered call
            s_dat[gt_col] = new_gt
            new_vcf.append(':'.join(s_dat))
    return new_vcf


def check_missing(vcf_line, max_miss):
    """Return boolean (True/False) for whether the VCF line passes filtering
    criteria based on missing data. A True means we keep the site and a False
    means we discard the site."""
    # Calculate the number of samples
    n_samp = len(vcf_line[9:])
    fmt_spec = vcf_line[8].split(':')
    gt_col = fmt_spec.index('GT')
    # Then count up the missing calls
    n_miss = 0
    for samp in vcf_line[9:]:
        if samp == '.':
            n_miss += 1
        else:
            s_dat = samp.split(':')
            if s_dat[gt_col] == '.':
                n_miss += 1
    # Calculate the proportion of missing data
    p_miss = float(n_miss) / float(n_samp)
    if p_miss > max_miss:
        return False
    return True


def main(vcf):
    """Main function"""
    with gzip.open(vcf, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                print(line.strip())
            else:
                vcf_rec = line.strip().split('\t')
                # First check that the REF and ALT fields are single characters
                ref = vcf_rec[3]
                alt = vcf_rec[4]
                # Skip the site if they are not - these are indels, "complex"
                # variants, or multiallelic sites.
                if len(ref) != 1 or len(alt) != 1:
                    continue
                # Next check the variant quality score. Skip the sites that are
                # below our threshold
                qual = float(vcf_rec[5])
                if qual < MIN_VAR_QUAL:
                    continue
                # Then, replace the low-depth calls with missing calls
                new_rec = geno_depth_flt(vcf_rec, MIN_GENO_DP)
                # Then, check the missing data fraction in this filtered
                # record
                keep = check_missing(new_rec, MAX_SITE_MISSING)
                if keep:
                    print('\t'.join(new_rec))
    return


main(vcf_in)
