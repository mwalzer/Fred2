# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: Variant
   :synopsis: Contains relevant classes describing variants.
.. moduleauthor:: schubert, walzer

"""
import math
from functools import total_ordering
from Fred2.Core.Base import MetadataLogger


# IMPORTANT order for integration precedence (due to genome position sorting and implementation in var.__lq__, smaller
# is of higher precedence): DEL before INS ... before SNP
VariationType = (lambda **enums: type('Enum', (), enums))(DEL=0,
                                                          INS=1,
                                                          FSDEL=2,
                                                          FSINS=3,
                                                          SNP=4,
                                                          UNKNOWN=5)
"""
Enum for variation types:
type.SNP, type.DEL, type.INS, type.FSDEL, type.FSINS, type.UNKNOWN
"""

@total_ordering
class MutationSyntax():
    """
    This class represents the mutation syntax of a variant and stores its 
    transcript and protein position

    :param str transID: The :class:`~Fred2.Core.Transcript.Transcript` id
    :param int transPos: The position of the :class:`~Fred2.Core.Variant.Variant` within the
                         :class:`~Fred2.Core.Transcript.Transcript`
    :param int protPos: The :class:`~Fred2.Core.Protein.Protein` position of the :class:`~Fred2.Core.Variant.Variant`
                        within the :class:`~Fred2.Core.Transcript.Transcript`
    :param str cds: The complete cds_mutation_syntax string
    :param str aas: The complete protein_mutation_syntax string
    """
    def __init__(self, transID, transPos, protPos, cds, aas, geneID=None):
        #TODO: is protPos always given? what about synonymous variants?
        self.transID = transID
        self.geneID = geneID
        self.tranPos = transPos
        self.protPos = protPos
        self.cdsMutationSyntax = cds  #c. ...
        self.aaMutationSyntax = aas  #p. ...

    def __eq__(self, other):
        return (self.cdsMutationSyntax, self.aaMutationSyntax, self.transID) == \
               (other.cdsMutationSyntax, other.aaMutationSyntax, other.transID)

    def __lt__(self, other):
        return (self.transID, self.tranPos) <\
               (other.transID, other.tranPos)


@total_ordering
class Variant(MetadataLogger):
    """
    A :class:`~Fred2.Core.Variant.Variant` contains information about a single genetic modification of
    the reference genome.

    """
    def __init__(self, id, type, chrom, genomePos, ref, obs, coding,
                 isHomozygous, isSynonymous, experimentalDesign=None, metadata=None):
        """
        Constructor for a variant, see init-types in class parameters

        :param str id: :class:`~Fred2.Core.Variant.Variant` id
        :param type: An Enum type of the :class:`~Fred2.Core.Variant.Variant` either SNP, DEL, or INS
        :type type: :func:`~Fred2.Core.Variant.VariationType`
        :param str chrom: The chromosome on which the variant lies
        :param int genomePos: The genomic position of the :class:`~Fred2.Core.Variant.Variant`
        :param str ref: The reference seq at the genomic position
        :param str obs: The observed variation at the genomic position
        :param coding: A dictionary of associated transcripts. Key=transcript_id,
                       value=:class:`~Fred2.Core.Variant.MutationSyntax`
        :type coding: dict(str,:class:`~Fred2.Core.Variant.MutationSyntax`)
        :param bool isHomozygous: Defines if variant is homozygous or not
        :param bool isSynonymous: Defines if variant is a synonymous mutation or not
        :param str experimentalDesign: String specifying the experimental condition (e.g. tumor)
        :param dict(list) metadata: meta information (not relevant for core functionality of Fred2)
        """
        MetadataLogger.__init__(self)
        self.id = id
        self.type = type
        self.chrom = chrom
        self.genomePos = genomePos
        self.ref = ref.upper()
        self.obs = obs.upper()
        self.isHomozygous = isHomozygous
        self.isSynonymous = isSynonymous
        self.coding = coding  # dict transcript_id:MutationSyntax
        self.experimentalDesign = "" if experimentalDesign is None else experimentalDesign

        if metadata is not None:
            for meta in metadata:
                self.log_metadata(meta, metadata[meta])

    def __repr__(self):
        return "Variant(g.%i%s>%s):%s" % (self.genomePos, self.ref, self.obs, self.experimentalDesign) \
            if self.experimentalDesign else "Variant(g.%i%s>%s)" % (self.genomePos, self.ref, self.obs)

    def __eq__(self, other):
        return self.id == other.id and self.coding == other.coding

    def __lt__(self, other):
        return (self.chrom, self.genomePos, self.type) <\
               (other.chrom, other.genomePos, other.type)

    def get_transcript_offset(self):
        """
        Returns the sequence offset caused by the mutation

        :return: The sequence offset
        :rtype: int
        """
        return len(self.obs) - len(self.ref)

    def get_shift(self):
        """
        Returns the frameshift offset caused by the mutation in {0,1,2}

        :return: The frameshift caused by mutation
        :rtype: int
        """
        return self.get_transcript_offset() % 3

    def get_annotated_transcript_pos(self, transID):
        """
        Returns the annotated :class:`~Fred2.Core.Transcript.Transcript` position

        :param str transID: The :class:`~Fred2.Core.Transcript.Transcript` ID of interest
        :return: The annotated :class:`~Fred2.Core.Transcript.Transcript` position of the given
                 :class:`~Fred2.Core.Transcript.Transcript` ID
        :rtype: int
        :raises KeyError: If variant is not annotated to the given :class:`~Fred2.Core.Transcript.Transcript` ID
        """
        trID = transID.split(":FRED_")[0]
        try:
            return self.coding[trID].tranPos
        except KeyError:
            raise KeyError("Variant {var} was not annotated to Transcript {tID}".format(var=repr(self), tID=transID))

    def get_annotated_protein_pos(self, transID):
        """
        Returns the annotated protein position

        :param str transID: The :class:`~Fred2.Core.Transcript.Transcript` ID of interest
        :return: The annotated :class:`~Fred2.Core.Protein.Protein` position of the given
                 :class:`~Fred2.Core.Transcript.Transcript` ID
        :rtype: int
        :raise KeyError: If :class:`~Fred2.Core.Variant.Variant` is not annotated to the given
                         :class:`~Fred2.Core.Transcript.Transcript` ID
        """
        trID = transID.split(":FRED_")[0]
        try:
            return self.coding[trID].protPos
        except KeyError:
            raise KeyError("Variant {var} was not annotated to " \
                           "Protein with transcript ID {tID}".format(var=repr(self),tID=transID))