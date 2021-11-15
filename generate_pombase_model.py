# from generate_rdf import GoCamModel, Annoton
# from gocamgen.gocamgen import GoCamModel, Annoton, GoCamEvidence
from ontobio.rdfgen.gocamgen.gocamgen import GoCamModel, Annoton, GoCamEvidence
from gaf_query import genes_and_annots_for_bp, NoCacheEagerRemoteSparqlOntology
from ontobio.ontol_factory import OntologyFactory
from ontobio.model.association import ymd_str
from pombase_golr_query import GeneConnectionSet, convert_relation_labels_to_curies, regulation_relations_curie, WITH_SUPPORT_FROM, AnnotationDataExtracter
from rdflib.term import URIRef
from ontobio.vocabulary.relations import OboRO
from prefixcommons.curie_util import expand_uri
import argparse

ro = OboRO()

ENABLED_BY = URIRef(expand_uri(ro.enabled_by))
PART_OF = URIRef(expand_uri(ro.part_of))
OCCURS_IN = URIRef(expand_uri(ro.occurs_in))
# HAS_INPUT = "has_input"
HAS_INPUT = convert_relation_labels_to_curies(["has_input"])[0]
# HAS_DIRECT_INPUT = "has_direct_input"
HAS_DIRECT_INPUT = convert_relation_labels_to_curies(["has_direct_input"])[0]
# DIRECTLY_REGULATES = "directly_regulates"
DIRECTLY_REGULATES = convert_relation_labels_to_curies(["directly_regulates"])[0]

# This is a hack. Prob don't need this after updating ontobio/../gocamgen:add_connection to handle Curie
RELATIONS_DICT = {
    "RO:0002400": "RO:0002400",
    "RO:0002233": "RO:0002233",
    # "has_regulation_target": "RO:0002211",  # regulates  # has_regulation_target is no longer in RO
    "RO:0011002": "RO:0002578",  # directly regulates
    "RO:0002614": "RO:0002233",  # has input
    "RO:0002578": "RO:0002578",
    "RO:0002629": "RO:0002629",
    "RO:0002630": "RO:0002630",
    "RO:0002325": "RO:0002325",
    "RO:0002326": "RO:0002326",
    "BFO:0000050": "BFO:0000050",
    "RO:0002263": "RO:0002263",
    "RO:0004035": "RO:0004035",
    "RO:0002264": "RO:0002264",
    "RO:0004034": "RO:0004034",
    "RO:0004033": "RO:0004033",
    "RO:0004032": "RO:0004032",
    "RO:0001025": "RO:0001025",
    "RO:0002432": "RO:0002432"  # is active in
}


class PomBaseAnnoton(Annoton):
    def __init__(self, gene_info, subject_id):
        self.enabled_by = subject_id
        self.molecular_function = self.get_aspect_object(gene_info, "molecular_function")
        self.cellular_component = self.get_aspect_object(gene_info, "cellular_component")
        self.biological_process = self.get_aspect_object(gene_info, "bp")
        self.connections = gene_info["connections"]
        self.individuals = {}

    #TODO: Should be replaced with gene_info.get(key)
    def get_aspect_object(self, gene_info, aspect):
        if aspect in gene_info:
            return gene_info[aspect]

parser = argparse.ArgumentParser()
parser.add_argument('-t', "--bp_term", type=str, required=True,
                    help="Biological process GO term that GOCAM should model. Comma-delimited list will create multiple model files")
parser.add_argument('-f', "--filename", type=str, required=False,
                    help="Destination filename - will end in '.ttl'")
parser.add_argument('-d', "--directory", type=str, required=False,
                    help="Destination directory - in this case titles will be auto-generated")
parser.add_argument('-g', "--gaf_source", type=str, required=True,
                    help="filename of GAF file to use as annotation source")
parser.add_argument('-o', "--go_ontology", type=str, required=False,
                    help="filename of GO ontology")
parser.add_argument('-j', "--tad_json", type=str, required=False,
                    help="Existing json data file to load into term-to-gene dictionary. Speeds up performance.")


def generate_model(bp_term, extracter: AnnotationDataExtracter, filename=None, directory=None, go_ontology=None):

    bp = bp_term
    if go_ontology is None:
        go_ontology = NoCacheEagerRemoteSparqlOntology("go")
    if not go_ontology.has_node(bp):
        raise "Invalid BP term {}".format(bp)

    gene_info = extracter.genes_and_annots_for_bp(bp_term)
    if gene_info is None:
        print("No gene_info for term:", bp_term)
        return

    model_title = "PomBase - " + bp + " - " + go_ontology.label(bp)
    if filename is None or directory is not None:
        # filename = model_title.replace(" - ", "_").replace(":", "_").replace(" ", "_")
        filename = model_title
        for unwanted_char in [" - ", ":", " ", "/", ","]:
            filename = filename.replace(unwanted_char, "_")

    model = GoCamModel(model_title, connection_relations=RELATIONS_DICT)
    model.writer.bp_id = model.declare_individual(bp)

    annotons = []
    for gene in gene_info:
        annoton = PomBaseAnnoton(gene_info[gene], gene)
        if annoton.molecular_function is not None:
            annotons.append(annoton)

    global_individuals_list = {}

    # translate lists of annotations
    for annoton in annotons:
        # writer.translate_annoton(annoton)
        sub = annoton.enabled_by
        obj = annoton.molecular_function.object

        # E.g. instance of gene product class. Keeping individuals at annoton level for our Pombase purposes.
        if sub not in annoton.individuals:
            enabler_id = model.declare_individual(sub)
            annoton.individuals[sub] = enabler_id
        else:
            enabler_id = annoton.individuals[sub]

        # E.g. instance of GO class
        object_id = str(obj.id)
        if object_id not in annoton.individuals:
            tgt_id = model.declare_individual(object_id)
            annoton.individuals[object_id] = tgt_id
        else:
            tgt_id = annoton.individuals[object_id]

        enabled_by_stmt = model.writer.emit(tgt_id, ENABLED_BY, enabler_id)
        part_of_stmt = model.writer.emit(tgt_id, PART_OF, model.writer.bp_id)

        axiom_id = model.add_axiom(enabled_by_stmt)
        association = annoton.molecular_function
        references = [str(ref) for ref in association.evidence.has_supporting_reference]
        evidence = GoCamEvidence(code=str(association.evidence.type),
                                 date=ymd_str(association.date, separator="-"),
                                 references=references)
        model.add_evidence(axiom_id, evidence)
        # Record annotation

        if annoton.cellular_component is not None:
            for cellular_component in annoton.cellular_component:
                cc_object_id = str(cellular_component.object.id)
                if cc_object_id not in annoton.individuals:
                    cc_id = model.declare_individual(cc_object_id)
                    annoton.individuals[cc_object_id] = cc_id
                else:
                    cc_id = annoton.individuals[cc_object_id]
                occurs_in_stmt = model.writer.emit(tgt_id, OCCURS_IN, cc_id)
                cc_axiom_id = model.add_axiom(occurs_in_stmt)
                cc_association = cellular_component
                cc_refs = [str(ref) for ref in cc_association.evidence.has_supporting_reference]
                cc_evidence = GoCamEvidence(code=str(cc_association.evidence.type),
                                            date=ymd_str(cc_association.date, separator="-"),
                                            references=cc_refs)
                model.add_evidence(cc_axiom_id, cc_evidence)
                # Record annotation
        
        global_individuals_list = {**global_individuals_list, **annoton.individuals}

    ### Connections - Now that all individuals should have been created
    global_connections_list = GeneConnectionSet()   # Tracking what connections have been added so far
    for annoton in annotons:
        for connection in annoton.connections.gene_connections:
            if connection.relation in [HAS_INPUT]:
                model.add_connection(connection, annoton)
                # Record annotation
                global_connections_list.append(connection)
            elif connection.relation in regulation_relations_curie + [HAS_DIRECT_INPUT]:
                if connection.relation == HAS_DIRECT_INPUT:
                    connection.relation = DIRECTLY_REGULATES
                # source_id = annoton.individuals[connection.gp_a]
                connection_gp_a = str(connection.gp_a)
                connection_gp_b = str(connection.gp_b)
                try:
                    source_id = annoton.individuals[str(connection.object_id)]
                except KeyError:
                    print("Say hey")
                    # New activity to declare for regulating annoton - also setup enabled_by triple and set source_id from MF
                    source_id = model.declare_individual(str(connection.object_id))
                    model.writer.emit(source_id, ENABLED_BY, annoton.individuals[connection_gp_a])
                    model.add_axiom((source_id, ENABLED_BY, annoton.individuals[connection_gp_a]))
                    annoton.individuals[str(connection.object_id)] = source_id
                # Probably need to switch source to be object (GO MF) of connection
                property_id = URIRef(expand_uri(model.connection_relations[str(connection.relation)]))
                # find annoton(s) of regulation target gene product
                target_annotons = model.writer.find_annotons(connection_gp_b, annotons)
                for t_annoton in target_annotons:
                    mf_annotation = t_annoton.get_aspect_object(gene_info[connection_gp_b], "molecular_function")
                    if mf_annotation is not None:
                        # target_id = global_individuals_list[mf_annotation["object"]["id"]]
                        target_id = t_annoton.individuals[str(mf_annotation.object.id)]
                        # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
                        model.writer.emit(source_id, property_id, target_id)
                        # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
                        model.writer.emit_axiom(source_id, property_id, target_id)
                        # Record mf_annotation
                global_connections_list.append(connection)

    # Now see if the with connections can fill anything in
    for annoton in annotons:
        for connection in annoton.connections.gene_connections:
            if connection.relation in [WITH_SUPPORT_FROM] and not global_connections_list.contains(connection):
                model.add_connection(connection, annoton)
                # Record annotation

    if directory is not None:
        filename = directory + filename
    model.write(filename)

    return model


if __name__ == "__main__":
    args = parser.parse_args()

    go_ontology = None
    if args.go_ontology:
        go_ontology = OntologyFactory().create(args.go_ontology)

    tad_json = None
    if args.tad_json:
        tad_json = args.tad_json

    extractor = AnnotationDataExtracter.parse_gaf(args.gaf_source, go_ontology, tad_json)

    for bp_term in args.bp_term.split(","):
        generate_model(bp_term, extractor, args.filename, args.directory, go_ontology=go_ontology)
