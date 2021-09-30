from ontobio.assoc_factory import AssociationSetFactory
from ontobio.io.gafparser import GafParser
from ontobio.io.assocparser import AssocParserConfig
from ontobio import ecomap


ecomapping = ecomap.EcoMap()


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
    acceptable_evidence_eco_codes = [ecomapping.coderef_to_ecoclass(ev) for ev in acceptable_evidence_codes]

    def filter_evidence(self, evidence_eco_codes=None):
        if evidence_eco_codes is None:
            evidence_eco_codes = self.acceptable_evidence_eco_codes
        new_gaf_list = []
        for a in self.gafs:
            if str(a.evidence.type) in evidence_eco_codes:
                new_gaf_list.append(a)
        self.gafs = new_gaf_list

    def __init__(self, filename, ontology, filter_evidence=True):
        parser_config = AssocParserConfig(ontology=ontology)
        parser = GafParser(config=parser_config)
        self.gafs = parser.parse(filename, skipheader=True)
        if filter_evidence:
            self.filter_evidence()
        self.association_set = AssociationSetFactory().create_from_assocs(self.gafs, ontology=ontology)

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
            if thing_key == "subject":
                thing = a.subject
            else:
                thing = a.object
            if str(thing.id) == thing_id and str(a.evidence.type) in GafAnnotationSet.acceptable_evidence_eco_codes:
                annots.append(a)
        return annots

    def subject_object_query(self, subject_id, object_id, gafs=None):
        return self.annotations_for_object(object_id, self.annotations_for_subject(subject_id))