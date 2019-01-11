from ontobio.assoc_factory import AssociationSetFactory
from ontobio.io.gafparser import GafParser

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

    def filter_evidence(self, evidence_codes=None):
        if evidence_codes is None:
            evidence_codes = GafAnnotationSet.acceptable_evidence_codes
        new_gaf_list = []
        for a in self.gafs:
            if a["evidence"]["type"] in evidence_codes:
                new_gaf_list.append(a)
        self.gafs = new_gaf_list

    def __init__(self, filename, ontology, filter_evidence=True):
        self.gafs = GafParser().parse(filename, skipheader=True)
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
            if a[thing_key]["id"] == thing_id and a["evidence"]["type"] in GafAnnotationSet.acceptable_evidence_codes:
                annots.append(a)
        return annots

    def subject_object_query(self, subject_id, object_id, gafs=None):
        return self.annotations_for_object(object_id, self.annotations_for_subject(subject_id))