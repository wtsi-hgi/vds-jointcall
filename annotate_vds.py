import argparse
import os

import hail as hl
import pyspark
from gnomad.utils.annotations import (
    annotate_adj,
    bi_allelic_site_inbreeding_expr,
)
from gnomad.utils.sparse_mt import (
    get_as_info_expr,
    get_site_info_expr,
)


def vds_to_annotated_mt(vds_path, mt_path, reference):
    vds = hl.vds.read_vds(vds_path)
    #filter to sites with 10 of fewer alleles.
    vd_rows = vds.variant_data.rows()
    sites_to_keep = vd_rows.filter(hl.len(vd_rows.alleles) <= 10)
    vds = hl.vds.filter_variants(vds, sites_to_keep)

    #make dense mt directly from the vds, then localize LGT/LAD/LPL/LPGT to
    #global allele indices (GT, AD, PL, PGT) on the densified entries
    hl.get_reference('GRCh38').add_sequence(reference)
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.annotate_entries(
        GT = hl.vds.lgt_to_gt(mt.LGT, mt.LA),
        AD = hl.vds.local_to_global(mt.LAD, mt.LA, hl.len(mt.alleles), fill_value=0, number='R'),
        PL = hl.vds.local_to_global(mt.LPL, mt.LA, hl.len(mt.alleles), fill_value=hl.missing('int32'), number='G'),
        PGT = hl.vds.lgt_to_gt(mt.LPGT, mt.LA),
    )

    #create qc_stat annotation for each site, using the dense mt's GT/AD/PL/DP/GQ
    qc_stat_ht = generate_qc_annotations(
        mt,
    )

    #remove hom ref only site
    mtf = mt.filter_rows(
        hl.agg.any(mt.GT.is_non_ref())
    )

    # extract annotatuon from gvcf info
    mt = mtf.transmute_entries(**mtf.gvcf_info)
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))
    
    # make info annptation

    #annotate DP, ReadPosRankSum, MQRankSum, BaseQRankSum, SB, MQ, FS, SOR
    #FS, SOR are bonus of SB annotation
    site_info_agg_fields = {
        'int32_sum_agg_fields': ['DP'],
        'median_agg_fields': ['ReadPosRankSum', 'MQRankSum', 'BaseQRankSum'],
        'array_sum_agg_fields': ['SB', 'RAW_MQandDP']
    }

    info_expr = get_site_info_expr(
        mt,
        sum_agg_fields=[],
        int32_sum_agg_fields=site_info_agg_fields["int32_sum_agg_fields"],
        median_agg_fields=site_info_agg_fields["median_agg_fields"],
        array_sum_agg_fields=site_info_agg_fields["array_sum_agg_fields"]
    )

    #annotate AS_MQ, AS_ReadPosRankSum, AS_MQRankSum, AS_QD, AS_SB_TABLE, AS_FS, AS_SOR
    # AS_FS, AS_SOR are bonus of AS_SB_TABLE annotation
    as_info_agg_fields = {
        'sum_agg_fields': ['AS_RAW_MQ'],
        'int32_sum_agg_fields': [],
        'median_agg_fields': ['AS_RAW_ReadPosRankSum', 'AS_RAW_MQRankSum', 'AS_QD'],
        'array_sum_agg_fields': ['AS_SB_TABLE'],
    }

    info_expr = info_expr.annotate(
        **get_as_info_expr(
            mt,
            sum_agg_fields=as_info_agg_fields["sum_agg_fields"],
            int32_sum_agg_fields=as_info_agg_fields["int32_sum_agg_fields"],
            median_agg_fields=as_info_agg_fields["median_agg_fields"],
            array_sum_agg_fields=as_info_agg_fields["array_sum_agg_fields"],
            treat_fields_as_allele_specific=True,
        )
    )

    # calculate DB:
    info_expr = info_expr.annotate(
        DB=hl.agg.any(hl.coalesce(mt.DB, False))
    )

    # calculate AC / AN / AF
    call_stats_expr = hl.agg.call_stats(mt.GT, mt.alleles)
    info_expr = info_expr.annotate(
        AN=call_stats_expr.AN,
        AC=call_stats_expr.AC[1:],
        AF=call_stats_expr.AF[1:],
    )

    info_ht = mt.select_rows(info=info_expr).rows()
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)


    # replace NaN with Na in FS and AS_FS
    mt = mt.annotate_rows(
        info = mt.info.annotate(
            FS=hl.or_missing(~hl.is_nan(mt.info.FS), mt.info.FS),
            AS_FS=mt.info.AS_FS.map(lambda x: hl.or_missing(~hl.is_nan(x), x)),
        )
    )

    #annotate QD
    mt = mt.annotate_rows(
        qual_stats=qc_stat_ht[mt.row_key].qc_stats
    )
    mt = mt.annotate_rows(
        info = mt.info.annotate(QD = mt.qual_stats.qd)
    )
    mt = mt.drop('qual_stats')
    
    #drop 'SB' and 'AS_SB_TABLE'
    mt = mt.annotate_rows(
        info = mt.info.drop('SB', 'AS_SB_TABLE')
    )
    #select only wanted entries
    mt = mt.select_entries(
        'AD',
        'DP',
        'GQ',
        'GT',
        'MIN_DP',
        'PGT',
        'PID',
        'PL',
        'PS',
        'RGQ',
        'SB'
    )
    # save mt to file
    print(f"Saving MT to path: {mt_path}")
    mt.write(mt_path, overwrite=True)


def generate_qc_annotations(mt: hl.MatrixTable) -> hl.Table:
    """
    Writes a Table with a per-site qd annotation, computed across all samples.

    The original version of this function also computed pab, gq, dp, nrq, ab,
    best_ab, nrdp, qual, combined_pab, and medians of several of those -- all
    removed here since only qd is actually used downstream (mt.info.QD).

    :param MatrixTable mt: Full MT
    :return: Table with qd annotation
    :rtype: Table
    """
    mt = mt.select_rows()

    lin = 10 ** (-mt.PL / 10.0)
    gp0 = lin[0] / hl.sum(lin)
    non_refs = mt.GT.is_non_ref()

    qual_agg = -10 * hl.agg.filter(
        non_refs,
        hl.agg.sum(hl.if_else(mt.PL[0] > 3000, -300, hl.log10(gp0)))
    )

    stats_expression = hl.struct(
        qd=hl.min(35, qual_agg / hl.agg.filter(non_refs, hl.agg.sum(mt.DP))),
    )

    return mt.annotate_rows(qc_stats=stats_expression).rows()


def file_uri(path):
    """Return path as an absolute file:// URI, stripping any existing file:// prefix."""
    if path.startswith("file://"):
        path = path[len("file://"):]
    return "file://" + os.path.abspath(path)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Densify a VDS, annotate site/AS info fields and write an annotated MatrixTable."
    )
    parser.add_argument("tmp_dir", help="Hail temporary directory")
    parser.add_argument("vds_path", help="Input VDS path")
    parser.add_argument("mt_path", help="Output MatrixTable path")
    parser.add_argument("reference", help="Reference FASTA path (GRCh38)")
    return parser.parse_args()


def main():
    args = parse_args()
    tmp_dir = file_uri(args.tmp_dir)
    vds_path = file_uri(args.vds_path)
    mt_path = file_uri(args.mt_path)
    reference = file_uri(args.reference)

    #initialise hail
    sc = pyspark.SparkContext()
    hadoop_config = sc._jsc.hadoopConfiguration()
    hl.init(sc=sc, tmp_dir=tmp_dir, default_reference="GRCh38")

    vds_to_annotated_mt(vds_path, mt_path, reference)

if __name__ == '__main__':
    main()