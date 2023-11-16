#!/usr/bin/env python3

import argparse
import csv
import json
import dataclasses
from typing import Dict, List, Set
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
parser.add_argument('-n', '--print_json', action='store_const', const=True)


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


@dataclasses.dataclass
class ModuleNode:
    ptn_id: str
    label: str
    terms: Set[str] = dataclasses.field(default_factory=set)
    leaf_genes: Set[str] = dataclasses.field(default_factory=set)

    def add_term(self, m_term: str):
        self.terms.add(m_term)

    def to_dict(self):
        return dataclasses.asdict(self)


@dataclasses.dataclass
class ProcessModule:
    term: str
    nodes: List[ModuleNode] = dataclasses.field(default_factory=list)

    def add_node(self, m_node: ModuleNode):
        existing_node = self.get_node_by_ptn_id(m_node.ptn_id)
        if existing_node and existing_node.leaf_genes == m_node.leaf_genes:
            existing_node.terms.update(m_node.terms)
        else:
            self.nodes.append(m_node)

    def get_node_by_ptn_id(self, ptn_id: str):
        for n in self.nodes:
            if n.ptn_id == ptn_id:
                return n
        return None

    def to_dict(self):
        return dataclasses.asdict(self)


@dataclasses.dataclass
class ModuleCollection:
    modules: List[ProcessModule] = dataclasses.field(default_factory=list)

    def module_by_term(self, m_term: str):
        for m in self.modules:
            if m.term == m_term:
                return m
        return None

    def add_module(self, m_dict: ProcessModule):
        self.modules.append(m_dict)

    def add_node_to_module(self, m_term: str, m_node: ModuleNode):
        m = self.module_by_term(m_term)
        if m:
            m.add_node(m_node)
        else:
            raise Exception("No module found for term '{}'".format(m_term))

    def to_dict(self):
        return dataclasses.asdict(self)


class SetHackJSONEncoder(json.JSONEncoder):
    # Can't believe json doesn't already convert sets to lists but whatevs.
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return super().default(obj)


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

    def get_human_ibas_by_ptn_and_term(ptn: str, term: str):
        if ptn in human_ibas_by_ptn and t in human_ibas_by_ptn[ptn]:
            iba_tuples = human_ibas_by_ptn[ptn][t]
            return iba_tuples
        else:
            return set()

    def print_iba_tuples(iba_tuples):
        if iba_tuples:
            [print("\t".join(list(iba))) for iba in iba_tuples]
        else:
            print("No human IBAs")


    modules = ModuleCollection()
    with open(args.clusters_file) as cf:
        cluster_reader = csv.reader(cf, delimiter="\t")
        if args.clusters_genes_is_edited:
            for r in cluster_reader:
                if r[0] in bp_gene_sets:
                    term = r[0]  # Set this for the next lines
                    if not args.print_json:
                        print("\t".join(r))  # Just print each line back out since we're only adding lines
                    bp_module = ProcessModule(term=term)
                    modules.add_module(bp_module)
                elif term and r[0].startswith("PTN"):
                    ptn = r[0]
                    t = r[2]
                    sf_label = ptn_to_sf_lkp[ptn]
                    term_label = ontology.label(t)
                    if not args.print_json:
                        print("\t".join([ptn, sf_label, t, term_label]))
                    leaf_genes = set()
                    if human_ibas_by_ptn:
                        leaf_gene_tuples = get_human_ibas_by_ptn_and_term(ptn, t)
                        leaf_genes = set([g[0] for g in leaf_gene_tuples])
                        if not args.print_json:
                            print_iba_tuples(leaf_gene_tuples)
                    node = ModuleNode(ptn_id=ptn, label=sf_label, terms={t}, leaf_genes=leaf_genes)
                    modules.add_node_to_module(term, node)
                elif not args.print_json:
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
                            leaf_gene_tuples = get_human_ibas_by_ptn_and_term(ptn, t)
                            if not args.print_json:
                                print_iba_tuples(leaf_gene_tuples)

    if args.print_json:
        modules_dict = dataclasses.asdict(modules)
        print(json.dumps(modules_dict, cls=SetHackJSONEncoder, indent=4))
