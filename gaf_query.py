from ontobio.assoc_factory import AssociationSetFactory
from ontobio.ontol_factory import OntologyFactory
from ontobio.io.gafparser import GafParser
from ontobio.sparql.sparql_ontology import EagerRemoteSparqlOntology, LazyRemoteSparqlOntology
from ontobio.ontol import Ontology
from pombase_direct_bp_annots_query import TermAnnotationDictionary, GOTermAnalyzer, ProgressTracker
from pombase_golr_query import AnnotationDataExtracter


class GafAnnotationSet():
    acceptable_evidence_codes = [
        "EXP",
        "IDA",
        "IPI",
        "IMP",
        "IGI",
        "IEP",
        "ND"
    ]

    def __init__(self, gaf_list):
        self.gafs = gaf_list

    def annotations_for_subject(self, subject_id, gafs=None):
        return self.annotations_for_thing("subject", subject_id, gafs)

    def annotations_for_object(self, object_id, gafs=None):
        return self.annotations_for_thing("object", object_id, gafs)

    def annotations_for_thing(self, thing_key, thing_id, gafs=None):
        annots = []
        if gafs is None:
            gafs = self.gafs
        for a in gafs:
            # if a[thing_key]["id"] == thing_id:
            if a[thing_key]["id"] == thing_id and a["evidence"]["type"] in GafAnnotationSet.acceptable_evidence_codes:
                annots.append(a)
        return annots

    def subject_object_query(self, subject_id, object_id, gafs=None):
        return self.annotations_for_object(object_id, self.annotations_for_subject(subject_id))

class NoCacheEagerRemoteSparqlOntology(EagerRemoteSparqlOntology):
    # Override cuz I don't care about text definitions or synonyms right now
    def subontology(self, nodes=None, **args):
        return Ontology.subontology(self, nodes, **args)

def genes_and_annots_for_bp(bp_term, filename, json_file=None):
    # ontology = OntologyFactory().create("go")
    ontology = NoCacheEagerRemoteSparqlOntology("go")
    afactory = AssociationSetFactory()
    a_set = afactory.create(ontology, file=filename, fmt="gaf")

    gafs = GafParser().parse(filename, skipheader=True)
    gas = GafAnnotationSet(gafs)

    tad = TermAnnotationDictionary(ontology, a_set, json_file)
    analyzer = GOTermAnalyzer(ontology)
    extracter = AnnotationDataExtracter(analyzer)
    
    gene_info = {}
    progress = ProgressTracker(len(tad.bps[bp_term]), "get relevant annotations/connections for each gene")
    for g in tad.bps[bp_term]:
        ### Find annots for g subject where go_term is MF - then check for extension part_of "GO:0010971"
        gene_annots = gas.annotations_for_subject(g)
        mf_annots = []
        for annot in gene_annots:
            if analyzer.is_molecular_function(annot["object"]["id"]):
                mf_annots.append(annot)
        gene_info[g] = {}
        gene_info[g]["bp"] = bp_term
        mf = extracter.relevant_mf_annotation(mf_annots, bp_term, tad) # 2a,b,c
        if mf is not None:
            gene_info[g]["molecular_function"] = mf
        cc = extracter.relevant_cc_annotation(mf_annots) # 2a
        if cc is not None:
            # ok_to_print_results = True
            gene_info[g]["cellular_component"] = [cc]
        cc_annots = []
        if "cellular_component" not in gene_info[g]:       
        #     cc_annots = extracter.direct_cc_annotations(gene_annots)
        #     if len(cc_annots) > 0:
        #         gene_info[g]["cellular_component"] = cc_annots[0] # 2b
            all_cc_annotations = extracter.cc_annotations(gene_annots)
            if all_cc_annotations and len(all_cc_annotations) == 1:
                gene_info[g]["cellular_component"] = all_cc_annotations
        gene_info[g]["connections"] = extracter.get_gene_connections(mf_annots, bp_term, tad)

        progress.print_progress()

    # extracter.print_results(gene_info, tad)
    return gene_info