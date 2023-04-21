import argparse
from ontobio import ontol_factory
from typing import Set
import os


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--term_file', help="Term list, typically a GO-slim")
parser.add_argument('-q', '--query_term_file', help="Term list to iterate through, checking if it is descendant of "
                                                    "terms in --term_file")
parser.add_argument('-o', '--ontology_file')
parser.add_argument('-d', '--destination_folder')


IS_A = "subClassOf"
PART_OF = "BFO:0000050"
REGULATES = "RO:0002211"
POSITIVELY_REGULATES = "RO:0002213"
NEGATIVELY_REGULATES = "RO:0002212"


def parse_term_list_file(term_file):
    terms = set()
    with open(term_file) as tf:
        for l in tf.readlines():
            terms.add(l.rstrip())
    return terms


def write_term_set_to_file(dest_file_path: str, term_set: Set):
    with open(dest_file_path, "w+") as df:
        for t in term_set:
            term_label = ontology.label(t)
            df.write("{}\n".format("\t".join([t, term_label])))


if __name__ == "__main__":
    args = parser.parse_args()

    ontology = ontol_factory.OntologyFactory().create(args.ontology_file)

    go_terms = parse_term_list_file(args.term_file)
    query_go_terms = parse_term_list_file(args.query_term_file)
    all_t_descendants_query_terms_intersection = set()
    print("term_file count:", len(go_terms))
    print("query_term_file count:", len(query_go_terms))

    slim_term_modules = {}
    go_term_files_written = 0
    for t in go_terms:
        # t_descendants = set(ontology.descendants(t, relations=[IS_A, PART_OF]))
        t_descendants = set(ontology.descendants(t, relations=[IS_A, PART_OF, REGULATES, POSITIVELY_REGULATES, NEGATIVELY_REGULATES]))
        t_descendants_query_terms_intersection = query_go_terms & t_descendants
        if t_descendants_query_terms_intersection:
            slim_term_modules[t] = t_descendants_query_terms_intersection
            all_t_descendants_query_terms_intersection = all_t_descendants_query_terms_intersection | t_descendants_query_terms_intersection

    cleaned_slim_term_modules = {}
    for t, module_terms in slim_term_modules.items():
        cleaned_slim_term_modules[t] = set()
        go_term_descendants_in_slim = set(ontology.descendants(t)) & (go_terms - set(t))
        if go_term_descendants_in_slim:
            module_terms_in_slim_descendants = set()
            for dt in go_term_descendants_in_slim:
                module_terms_in_slim_descendants = module_terms_in_slim_descendants | slim_term_modules.get(dt, set())
            cleaned_slim_term_modules[t] = slim_term_modules[t] - module_terms_in_slim_descendants
        else:
            cleaned_slim_term_modules[t] = slim_term_modules[t]

    for ct, module_terms in cleaned_slim_term_modules.items():
        out_filename = os.path.join(args.destination_folder, "{}.tsv".format(ct.replace(":", "_")))
        write_term_set_to_file(out_filename, module_terms)
        go_term_files_written += 1
    print("go_term_files_written:", go_term_files_written)

    # Terms in go_terms not having any query_go_terms in its descendants
    terms_that_dont_intersect = query_go_terms - all_t_descendants_query_terms_intersection
    print("terms_that_dont_intersect count:", len(terms_that_dont_intersect))
    if terms_that_dont_intersect:
        out_filename = os.path.join(args.destination_folder, "BP_module_terms_dont_trace_to_slim.tsv")
        write_term_set_to_file(out_filename, terms_that_dont_intersect)
