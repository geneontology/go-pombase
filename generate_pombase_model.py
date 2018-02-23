from generate_rdf import GoCamModel, Annoton
from gaf_query import genes_and_annots_for_bp
from pombase_golr_query import GeneConnectionSet
from rdflib.term import URIRef
from ontobio.vocabulary.relations import OboRO
from prefixcommons.curie_util import expand_uri
import argparse

ro = OboRO()

ENABLED_BY = URIRef(expand_uri(ro.enabled_by))
PART_OF = URIRef(expand_uri(ro.part_of))
OCCURS_IN = URIRef(expand_uri(ro.occurs_in))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', "--bp_term", type=str, required=True,
                        help="Biological process GO term that GOCAM should model")
    parser.add_argument('-f', "--filename", type=str, required=False,
                        help="Destination filename - will end in '.ttl'")
    parser.add_argument('-g', "--gaf_source", type=str, required=True,
                        help="filename of GAF file to use as annotation source")
    parser.add_argument('-j', "--tad_json", type=str, required=False,
                        help="Existing json data file to load into term-to-gene dictionary. Speeds up performance.")

    args = parser.parse_args()

    bp = args.bp_term
    gene_info = genes_and_annots_for_bp(bp, args.gaf_source, json_file=args.tad_json)

    model = GoCamModel(args.filename)
    model.writer.bp_id = model.declare_individual(bp)

    annotons = []
    for gene in gene_info:
        annoton = Annoton(gene_info[gene], gene)
        if annoton.molecular_function is not None:
            annotons.append(annoton)

    global_individuals_list = {}

    # translate lists of annotations
    for annoton in annotons:
        # Class
        if annoton.enabled_by not in model.classes:
            model.declare_class(annoton.enabled_by)

        # writer.translate_annoton(annoton)
        sub = annoton.enabled_by
        obj = annoton.molecular_function["object"]

        # E.g. instance of gene product class. Keeping individuals at annoton level for our Pombase purposes.
        if sub not in annoton.individuals:
            enabler_id = model.declare_individual(sub)
            annoton.individuals[sub] = enabler_id
        else:
            enabler_id = annoton.individuals[sub]

        # E.g. instance of GO class
        if obj["id"] not in annoton.individuals:
            tgt_id = model.declare_individual(obj["id"])
            annoton.individuals[obj["id"]] = tgt_id
        else:
            tgt_id = annoton.individuals[obj["id"]]

        enabled_by_stmt = model.writer.emit(tgt_id, ENABLED_BY, enabler_id)
        part_of_stmt = model.writer.emit(tgt_id, PART_OF, model.writer.bp_id)

        axiom_id = model.add_axiom(enabled_by_stmt)
        association = annoton.molecular_function
        model.add_evidence(axiom_id, association["evidence"]["type"], association["evidence"]["has_supporting_reference"])

        if annoton.cellular_component is not None:
            cc_object_id = annoton.cellular_component["object"]["id"]
            if cc_object_id not in annoton.individuals:
                cc_id = model.declare_individual(cc_object_id)
                annoton.individuals[cc_object_id] = cc_id
            occurs_in_stmt = model.writer.emit(tgt_id, OCCURS_IN, cc_id)
            cc_axiom_id = model.add_axiom(occurs_in_stmt)
            cc_association = annoton.cellular_component
            model.add_evidence(cc_axiom_id, cc_association["evidence"]["type"], cc_association["evidence"]["has_supporting_reference"])
        
        global_individuals_list = {**global_individuals_list, **annoton.individuals}

    ### Connections - Now that all individuals should have been created
    connection_relations = {"has_direct_input" : "RO:0002400",
                            "has input" : "RO:0002233",
                            "has_regulation_target" : "RO:0002211", # regulates
                            "regulates_activity_of" : "RO:0002578", # directly regulates
                            "with_support_from" : "RO:0002233" # has input
                            }
    global_connections_list = GeneConnectionSet()   # Tracking what connections have been added so far
    for annoton in annotons:
        for connection in annoton.connections.gene_connections:
            if connection.relation in ["has_direct_input", "has input"]:
                try:
                    source_id = annoton.individuals[connection.object_id]
                except KeyError:
                    source_id = model.declare_individual(connection.object_id)
                    annoton.individuals[connection.object_id] = source_id
                # object_source_id = annoton.individuals[annoton.molecular_function["object"]["id"]]
                # Add enabled by stmt for object_id - this is essentially adding another annoton connecting gene-to-extension/with-MF to the model
                model.writer.emit(source_id, ENABLED_BY, annoton.individuals[annoton.enabled_by])
                model.writer.emit_axiom(source_id, ENABLED_BY, annoton.individuals[annoton.enabled_by])
                property_id = URIRef(expand_uri(connection_relations[connection.relation]))
                target_id = global_individuals_list[connection.gp_b]
                # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
                model.writer.emit(source_id, property_id, target_id)
                # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
                model.writer.emit_axiom(source_id, property_id, target_id)
                global_connections_list.append(connection)
            elif connection.relation in ["has_regulation_target", "regulates_activity_of"]:
                # source_id = annoton.individuals[connection.gp_a]
                source_id = annoton.individuals[connection.object_id]
                # Probably need to switch source to be object (GO MF) of connection
                property_id = URIRef(expand_uri(connection_relations[connection.relation]))
                # find annoton(s) of regulation target gene product
                target_annotons = model.writer.find_annotons(connection.gp_b, annotons)
                for t_annoton in target_annotons:
                    mf_annotation = t_annoton.get_aspect_object(gene_info[connection.gp_b], "molecular_function")
                    if mf_annotation is not None:
                        target_id = global_individuals_list[mf_annotation["object"]["id"]]
                        # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
                        model.writer.emit(source_id, property_id, target_id)
                        # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
                        model.writer.emit_axiom(source_id, property_id, target_id)
                global_connections_list.append(connection)

    # Now see if the with connections can fill anything in
    for annoton in annotons:
        for connection in annoton.connections.gene_connections:
            if connection.relation in ["with_support_from"] and not global_connections_list.contains(connection):
                try:
                    source_id = annoton.individuals[connection.object_id]
                except KeyError:
                    source_id = model.declare_individual(connection.object_id)
                    annoton.individuals[connection.object_id] = source_id
                model.writer.emit(source_id, ENABLED_BY, annoton.individuals[annoton.enabled_by])
                model.writer.emit_axiom(source_id, ENABLED_BY, annoton.individuals[annoton.enabled_by])
                property_id = URIRef(expand_uri(connection_relations[connection.relation]))
                target_id = global_individuals_list[connection.gp_b]
                # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
                model.writer.emit(source_id, property_id, target_id)
                # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
                model.writer.emit_axiom(source_id, property_id, target_id)

    with open(model.modeltitle + ".ttl", 'wb') as f:
        model.writer.writer.serialize(destination=f)

if __name__ == "__main__":
    main()