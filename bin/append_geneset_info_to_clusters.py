#!/usr/bin/env python3

import argparse
import csv
from typing import Dict
from go_bp_modeling.pombase_direct_bp_annots_query import TermAnnotationDictionary
from ontobio.ontol_factory import OntologyFactory


parser = argparse.ArgumentParser()
parser.add_argument('-j', '--tad_json')
parser.add_argument('-c', '--clusters_file')
parser.add_argument('-s', '--ptn_to_sf_file')
parser.add_argument('-o', '--ontology_file')
parser.add_argument('-g', '--iba_gaf_file')
parser.add_argument('-e', '--clusters_genes_is_edited', action='store_const', const=True,
                    help="If flag set, --clusters_file will be parsed at term and gene level (as opposed to just term)")


def parse_ptn_to_sf(ptn_to_sf_file):
    ptn_to_sf = {}
    with open(ptn_to_sf_file) as psf:
        reader = csv.reader(psf, delimiter="\t")
        for r in reader:
            an_id = r[0]
            ptn = r[1]
            sf_id = r[2]
            sf_label = r[3]
            ptn_to_sf[ptn] = sf_label
    return ptn_to_sf


if __name__ == "__main__":
    args = parser.parse_args()
    bp_gene_sets = TermAnnotationDictionary.parse_json_to_bp_gene_sets(args.tad_json)
    ptn_to_sf_lkp = parse_ptn_to_sf(args.ptn_to_sf_file)
    ontology = OntologyFactory().create(args.ontology_file)
    human_ibas_by_ptn = {}
    if args.iba_gaf_file:
        with open(args.iba_gaf_file) as gf:
            reader = csv.reader(gf, delimiter="\t")
            for r in reader:
                if r[0].startswith("!"):
                    # header line
                    continue
                gene_symbol = r[2]
                go_term = r[4]
                with_from_col = r[7]
                db_object_synonyms = r[10]

                pthr_ibd_ptn, the_rest = with_from_col.split("|", maxsplit=1)
                pthr, ibd_ptn = pthr_ibd_ptn.split(":", maxsplit=1)
                gene_uniprot_id, leaf_ptn = db_object_synonyms.split("|", maxsplit=1)
                iba_data_tuple = (gene_uniprot_id, gene_symbol)

                if ibd_ptn not in human_ibas_by_ptn:
                    human_ibas_by_ptn[ibd_ptn] = {}
                if go_term not in human_ibas_by_ptn[ibd_ptn]:
                    human_ibas_by_ptn[ibd_ptn][go_term] = set()
                human_ibas_by_ptn[ibd_ptn][go_term].add(iba_data_tuple)

    def print_human_ibas_by_ptn_and_term(ptn: str, term: str):
        if ptn in human_ibas_by_ptn and t in human_ibas_by_ptn[ptn]:
            iba_tuples = human_ibas_by_ptn[ptn][t]
            [print("\t".join(list(iba))) for iba in iba_tuples]
        else:
            print("No human IBAs")

    with open(args.clusters_file) as cf:
        cluster_reader = csv.reader(cf, delimiter="\t")
        if args.clusters_genes_is_edited:
            for r in cluster_reader:
                if r[0] in bp_gene_sets:
                    term = r[0]  # Set this for the next lines
                    print("\t".join(r))  # Just print each line back out since we're only adding lines
                # elif r[0].startswith("#"):
                #     print("\t".join(r))  # Just print comments back out
                elif term and r[0].startswith("PTN"):
                    # Hoping this line
                    ptn = r[0]
                    t = r[2]
                    # genes_and_terms = bp_gene_sets[term]
                    # if ("PANTHER:{}".format(ptn), t) in bp_gene_sets[term]:
                    sf_label = ptn_to_sf_lkp[ptn]
                    term_label = ontology.label(t)
                    print("\t".join([ptn, sf_label, t, term_label]))
                    if human_ibas_by_ptn:
                        print_human_ibas_by_ptn_and_term(ptn, t)
                    # else:
                    #     print("Nothing printed")
                else:
                    print("\t".join(r))  # Nothing cool. Just print the line back out.
        else:
            for r in cluster_reader:
                print("\t".join(r))  # Just print each line back out since we're only adding lines
                if r[0] in bp_gene_sets:
                    term = r[0]
                    genes_and_terms = bp_gene_sets[term]
                    for g, t in genes_and_terms:
                        ptn = g.replace("PANTHER:", "")
                        sf_label = ptn_to_sf_lkp[ptn]
                        term_label = ontology.label(t)
                        print("\t".join([ptn, sf_label, t, term_label]))
                        if human_ibas_by_ptn:
                            print_human_ibas_by_ptn_and_term(ptn, t)
