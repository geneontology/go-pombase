from ontobio.assoc_factory import AssociationSetFactory
from ontobio.ontol_factory import OntologyFactory
from ontobio.io.gafparser import GafParser
from ontobio.sparql.sparql_ontology import EagerRemoteSparqlOntology, LazyRemoteSparqlOntology
from ontobio.ontol import Ontology
from pombase_direct_bp_annots_query import TermAnnotationDictionary, GOTermAnalyzer, ProgressTracker
from pombase_golr_query import AnnotationDataExtracter
from gaf_annotation_set import GafAnnotationSet

class NoCacheEagerRemoteSparqlOntology(EagerRemoteSparqlOntology):
    # Override cuz I don't care about text definitions or synonyms right now
    def subontology(self, nodes=None, **args):
        return Ontology.subontology(self, nodes, **args)

def genes_and_annots_for_bp(bp_term, filename, json_file=None, go_ontology=None):
    if go_ontology is None:
        go_ontology = NoCacheEagerRemoteSparqlOntology("go")
    # afactory = AssociationSetFactory()
    # a_set = afactory.create(ontology, file=filename, fmt="gaf")

    # gafs = GafParser().parse(filename, skipheader=True)
    gas = GafAnnotationSet(filename, go_ontology, filter_evidence=True)
    gas.filter_evidence()

    tad = TermAnnotationDictionary(go_ontology, gas.association_set, json_file)
    analyzer = GOTermAnalyzer(go_ontology)
    extracter = AnnotationDataExtracter(analyzer)
    
    gene_info = {}
    progress = ProgressTracker(len(tad.bps[bp_term]), "get relevant annotations/connections for each gene")
    for g in tad.bps[bp_term]:
        ### Find annots for g subject where go_term is MF - then check for extension part_of "GO:0010971"
        gene_annots = gas.annotations_for_subject(g)
        mf_annots = []
        for annot in gene_annots:
            if analyzer.is_molecular_function(str(annot.object.id)):
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
        if "cellular_component" not in gene_info[g]:
            all_cc_annotations = extracter.cc_annotations(gene_annots)
            if all_cc_annotations and len(all_cc_annotations) == 1:
                gene_info[g]["cellular_component"] = all_cc_annotations
        gene_info[g]["connections"] = extracter.get_gene_connections(mf_annots, bp_term, tad)

        progress.print_progress()

    # extracter.print_results(gene_info, tad)
    return gene_info