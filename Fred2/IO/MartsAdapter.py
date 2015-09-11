# This code is part of the Fred2 distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
.. module:: IO.MartsAdapter
   :synopsis: BDB-Adapter for BioMart
.. moduleauthor:: walzer, schubert
"""

import csv
import urllib2
import warnings
import logging

from Fred2.IO.ADBAdapter import ADBAdapter, EAdapterFields
import MySQLdb

class MartsAdapter(ADBAdapter):
    def __init__(self, usr=None, host=None, pwd=None, db=None, biomart=None):
        """
        used to fetch sequences from given RefSeq id's either from BioMart if no credentials given else from a MySQLdb
        :param usr: db user e.g. = 'ucsc_annot_query'
        :param host: db host e.g. = "pride"
        :param pwd: pw for user e.g. = 'an0q3ry'
        :param db: db on host e.g. = "hg18_ucsc_annotation"
        """
        self.ids_proxy = dict()
        self.gene_proxy = dict()
        self.sequence_proxy = dict()

        if usr and host and pwd and db:
            self.connection = MySQLdb.connect(user=usr, host=host, passwd=pwd, db=db)
        else:
            self.connection = None

        if biomart:
            self.biomart_url = biomart
            if not self.biomart_url.endswith("/biomart/martservice?query="):
                self.biomart_url += "/biomart/martservice?query="
        else:
            self.biomart_url = "http://biomart.org/biomart/martservice?query="
            #new: http://central.biomart.org/biomart/martservice?query="
            #http://www.ensembl.org/biomart/martview/
            #http://grch37.ensembl.org/biomart/
            #http://localhost:9000/biomart/martservice?query=%3C!DOCTYPE%20Query%3E%3CQuery%20client=%22biomartclient%22%20processor=%22TSV%22%20limit=%22-1%22%20header=%221%22%3E%3CDataset%20name=%22hsapiens_gene_ensembl%22%20config=%22gene_ensembl_config%22%3E%3CFilter%20name=%22uniprot_genename%22%20value=%22TP53%22%20filter_list=%22%22/%3E%3CFilter%20name=%22germ_line_variation_source%22%20value=%22dbSNP%22%20filter_list=%22%22/%3E%3CAttribute%20name=%22snp_ensembl_gene_id%22/%3E%3CAttribute%20name=%22snp_chromosome_name%22/%3E%3CAttribute%20name=%22snp_ensembl_transcript_id%22/%3E%3CAttribute%20name=%22snp_start_position%22/%3E%3CAttribute%20name=%22snp_ensembl_peptide_id%22/%3E%3CAttribute%20name=%22snp_end_position%22/%3E%3CAttribute%20name=%22snp_external_gene_name%22/%3E%3CAttribute%20name=%22snp_strand%22/%3E%3CAttribute%20name=%22source_description%22/%3E%3CAttribute%20name=%22germ_line_variation_source%22/%3E%3CAttribute%20name=%22allele%22/%3E%3CAttribute%20name=%22variation_name%22/%3E%3C/Dataset%3E%3C/Query%3E
            #self.biomart_url = """http://134.2.9.124/biomart/martservice?query="""
        self.biomart_head = """
        <?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query client="true" processor="TSV" limit="-1" header="1" uniqueRows = "1" >
                <Dataset name="%s" config="%s">
        """.strip()
        self.biomart_tail = """
                </Dataset>
            </Query>
        """.strip()
        self.biomart_filter = """<Filter name="%s" value="%s" filter_list=""/>"""
        self.biomart_attribute = """<Attribute name="%s"/>"""

    #TODO: refactor ... function based on old code
    def get_product_sequence(self, product_refseq, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        fetches product sequence for the given id
        :param product_refseq: given refseq id
        :return: list of dictionaries of the requested sequence, the respective strand and the associated gene name
        """

        if product_refseq in self.sequence_proxy:
            return self.sequence_proxy[product_refseq]

        if self.connection:
            cur = self.connection.cursor()
            #query = "SELECT * FROM sbs_ncbi_mrna WHERE id='%s';"%('5')
            query = "SELECT refseq,seq FROM hg19_ucsc_annotation.refLink INNER JOIN sbs_ncbi_protein on protAcc=refseq WHERE mrnaAcc='%s';" % (
                product_refseq)  # gives NP plus prot seq
            #mysql -h genome-mysql.cse.ucsc.edu -u genome -D hg19 -N -A -e 'SELECT refseq,seq FROM hg19_ucsc_annotation.refLink INNER JOIN sbs_ncbi_protein on protAcc=refseq WHERE mrnaAcc='NM_002486';'
            try:
                cur.execute(query)
                result = cur.fetchall()
                if result:
                    product_refseq = result[0][0]
                    product_sequence = result[0][1]
                else:
                    warnings.warn("An Error occured while executing query: " + query + "\n")
                    return None
            except MySQLdb.Error, message:
                warnings.warn(
                    "An Error occured while executing query: " + query + "\n" + "Error message: " + message[1])
                return None
            self.sequence_proxy[product_refseq] = product_sequence
            return [{product_refseq: product_sequence}]
        else:
            filter = None
            if product_refseq.startswith('NP_'):
                filter = "refseq_peptide"
            elif product_refseq.startswith('XP_'):
                filter = "refseq_peptide_predicted"
            elif product_refseq.startswith('ENS'):
                filter = "ensembl_peptide_id"
            else:
                warnings.warn("No correct product id: " + product_refseq)
                return None
            rq_n = self.biomart_head%(_db,_dataset) \
                + self.biomart_filter%(filter, str(product_refseq))  \
                + self.biomart_attribute%("peptide")  \
                + self.biomart_attribute%(filter)  \
                + self.biomart_attribute%("external_gene_id")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_tail

            logging.warn(rq_n)

            tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
            tsvselect = [x for x in tsvreader]
            if not tsvselect:
                return None
            self.sequence_proxy[product_refseq] = tsvselect[0]["Protein"]
            return self.sequence_proxy[product_refseq]  # TODO to SeqRecord

    #TODO: refactor ... function based on old code
    def get_transcript_sequence(self, transcript_refseq, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        Fetches transcript sequence for the given id
        :param transcript_refseq:
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        # TODO transcript_refseq to transcript_id, sniff which, query according
        if transcript_refseq in self.sequence_proxy:
            return self.sequence_proxy[transcript_refseq]

        if self.connection:
            cur = self.connection.cursor()
            #query = "SELECT * FROM sbs_ncbi_mrna WHERE id='%s';"%('5')
            query = "SELECT refseq,seq FROM hg19_ucsc_annotation.sbs_ncbi_mrna WHERE refseq LIKE '" + transcript_refseq + "%'"
            try:
                cur.execute(query)
                result = cur.fetchall()
                if result:
                    transcript_refseq = result[0][0]
                    transcript_sequence = result[0][1]
                    if len(result) > 1:
                        warnings.warn("Ambiguous transcript refseq query: " + transcript_refseq + "\n")
                else:
                    warnings.warn("An Error occured while executing query: " + query + "\n")
                    return None
            except MySQLdb.Error, message:
                warnings.warn("An Error occured while executing query: " + query + "\n" + "Error message: " + message[1])
                return None
            self.sequence_proxy[transcript_refseq] = transcript_sequence
            return [{transcript_refseq: transcript_sequence}]
        else:
            filter = None
            if transcript_refseq.startswith('NM_'):
                filter = "refseq_mrna"
            elif transcript_refseq.startswith('XM_'):
                filter = "refseq_mrna_predicted"
            elif transcript_refseq.startswith('ENS'):
                filter = "ensembl_transcript_id"
            else:
                warnings.warn("No correct transcript id: " + transcript_refseq)
                return None
            rq_n = self.biomart_head%(_db, _dataset) \
                + self.biomart_filter%(filter, str(transcript_refseq))  \
                + self.biomart_attribute%(filter)  \
                + self.biomart_attribute%("coding")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_tail

            tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
            tsvselect = [x for x in tsvreader]
            self.sequence_proxy[transcript_refseq] = tsvselect[0]['Coding sequence']  # TODO is strand information worth to pass along?
            return self.sequence_proxy[transcript_refseq]

    def get_transcript_information(self, transcript_refseq, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        It also already uses the Field-Enum for DBAdapters

        Fetches transcript sequence for the given id
        :param transcript_refseq:
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        # TODO transcript_refseq to transcript_id, sniff which, query according
        if transcript_refseq in self.sequence_proxy:
            return self.sequence_proxy[transcript_refseq]

        if self.connection:
            cur = self.connection.cursor()
            #query = "SELECT * FROM sbs_ncbi_mrna WHERE id='%s';"%('5')
            query = "SELECT refseq,seq FROM hg19_ucsc_annotation.sbs_ncbi_mrna WHERE refseq LIKE '" + transcript_refseq + "%'"
            try:
                cur.execute(query)
                result = cur.fetchall()
                if result:
                    transcript_refseq = result[0][0]
                    transcript_sequence = result[0][1]
                    if len(result) > 1:
                        warnings.warn("Ambiguous transcript refseq query: " + transcript_refseq + "\n")
                else:
                    warnings.warn("An Error occured while executing query: " + query + "\n")
                    return None
            except MySQLdb.Error, message:
                warnings.warn("An Error occured while executing query: " + query + "\n" + "Error message: " + message[1])
                return None
            self.ids_proxy[transcript_refseq] = transcript_sequence
            return [{transcript_refseq: transcript_sequence}]
        else:
            filter = None
            if transcript_refseq.startswith('NM_'):
                filter = "refseq_mrna"
            elif transcript_refseq.startswith('XM_'):
                filter = "refseq_mrna_predicted"
            elif transcript_refseq.startswith('ENS'):
                filter = "ensembl_transcript_id"
            else:
                warnings.warn("No correct transcript id: " + transcript_refseq)
                return None
            rq_n = self.biomart_head%(_db, _dataset) \
                + self.biomart_filter%(filter, str(transcript_refseq))  \
                + self.biomart_attribute%(filter)  \
                + self.biomart_attribute%("coding")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_tail

            print "Transcript information ",rq_n
            tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
            tsvselect = [x for x in tsvreader]
            if not tsvselect:
                warnings.warn("No entry found for ID %s"%transcript_refseq)
                return None

            self.ids_proxy[transcript_refseq] = {EAdapterFields.SEQ: tsvselect[0]['Coding sequence'],
                                                      EAdapterFields.GENE: tsvselect[0].get('Associated Gene Name', ""),
                                                      EAdapterFields.STRAND: "-" if int(tsvselect[0]['Strand']) < 0
                                                      else "+"}
            return self.ids_proxy[transcript_refseq]

    def get_transcript_position(self, start, stop, gene_id, transcript_id, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        If no transcript position is available for the variant
        :param start:
        :param stop:
        :param gene_id:
        :param transcript_id:
        :param _db:
        :param _dataset:
        :return:
        """
        # ma = MartsAdapter(biomart="http://grch37.ensembl.org")
        # print ma.get_transcript_position('17953929', '17953943', 'ENSG00000074964', 'ENST00000361221')
        # (1674, 1688)
        if start + stop + gene_id + transcript_id in self.gene_proxy:
            return self.gene_proxy[start + stop + gene_id + transcript_id]

        try:
            x = int(start)
            y = int(stop)
        except Exception as e:
            logging.warning(','.join([str(start), str(stop)]) + ' does not seem to be a genomic position.')
            return None

            filter_g = None
            filter_t = None
            if transcript_id.startswith('NM_'):
                filter_t = "refseq_mrna"
            elif transcript_id.startswith('XM_'):
                filter_t = "refseq_mrna_predicted"
            elif transcript_id.startswith('ENS'):
                filter_t = "ensembl_transcript_id"
            else:
                warnings.warn("No correct element id: " + transcript_id)
                return None
            if gene_id.startswith('ENS'):
                filter_g = "ensembl_gene_id"
            else:
                filter_g = "uniprot_genename"
                warnings.warn("Could not determine the type og gene_id: " + gene_id)

        rq_n = self.biomart_head%(_db, _dataset) \
            + self.biomart_filter%("ensembl_gene_id", str(gene_id))  \
            + self.biomart_filter%("ensembl_transcript_id", str(transcript_id))  \
            + self.biomart_attribute%("exon_chrom_start")  \
            + self.biomart_attribute%("exon_chrom_end")  \
            + self.biomart_tail

        tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url + urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        exons = [ex for ex in tsvreader]

        pos_sum = 0
        if not exons:
            logging.warning(','.join([str(gene_id), str(transcript_id)]) + ' does not seem to have exons mapped.')
            return None
        for e in exons:
            se = int(e["Exon Chr Start (bp)"])
            ee = int(e["Exon Chr End (bp)"])
            print se, ee + 1
            if x in range(se, ee + 1):
                if not y in range(se, ee + 1):
                    logging.warning(','.join([str(start), str(stop)]) + ' seems to span more than one exon.')
                    return None
                else:
                    self.gene_proxy[start + stop + gene_id + transcript_id] = (x - se + 1 + pos_sum, y - se + 1 + pos_sum)
                    return x - se + 1 + pos_sum, y - se + 1 + pos_sum
            else:
                pos_sum += ee - se + 1



    #TODO: refactor ... function based on old code
    def get_variant_gene(self, chrom, start, stop, _db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        Fetches the important db ids and names for given chromosomal location
        :param chrom: integer value of the chromosome in question
        :param start: integer value of the variation start position on given chromosome
        :param stop: integer value of the variation stop position on given chromosome
        :return: The respective gene name, i.e. the first one reported
        """
        if chrom + start + stop in self.gene_proxy:
            return self.gene_proxy[chrom + start + stop]

        rq_n = self.biomart_head%(_db, _dataset) \
            + self.biomart_filter%("chromosome_name", str(chrom))  \
            + self.biomart_filter%("start", str(start))  \
            + self.biomart_filter%("end", str(stop))  \
            + self.biomart_attribute%("uniprot_genename")  \
            + self.biomart_tail

        tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if tsvselect and tsvselect[0]:
            self.gene_proxy[chrom + start + stop] = tsvselect[0]['UniProt Gene Name']
            return tsvselect[0]['UniProt Gene Name']
        else:
            logging.warning(','.join([str(chrom), str(start), str(stop)]) + ' does not denote a known gene location')
            return ''

    def get_transcript_information_from_protein_id(self, **kwargs):
        """
        It also already uses the Field-Enum for DBAdapters

        Fetches transcript sequence for the given id
        :param transcript_refseq:
        :return: list of dictionary of the requested sequence, the respective strand and the associated gene name
        """
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        filter = None
        db_id = ""
        if "refseq" in kwargs:
            filter = "refseq_peptide"
            db_id = kwargs["refseq"]
        elif "ensemble" in kwargs:
            filter = "ensembl_peptide_id"
            db_id = kwargs["ensemble"]
        elif "swiss_accid" in kwargs:
            filter = "uniprot_swissprot_accession"
            db_id = kwargs["swiss_accid"]
        elif "swiss_gene" in kwargs:
            filter= "uniprot_swissprot"
            db_id = kwargs["swiss_gene"]
        else:
            warnings.warn("No correct transcript id")
            return None
        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(filter, str(db_id)) \
               + self.biomart_attribute%(filter) \
               + self.biomart_attribute%("coding") \
               + self.biomart_attribute%("external_gene_id") \
               + self.biomart_attribute%("strand") \
               + self.biomart_tail

        tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if not tsvselect:
            warnings.warn("No entry found for ID %s"%db_id)
            return None

        self.ids_proxy[db_id] = {EAdapterFields.SEQ: tsvselect[0]['Coding sequence'],
                                             EAdapterFields.GENE: tsvselect[0]['Associated Gene Name'],
                                             EAdapterFields.STRAND: "-" if int(tsvselect[0]['Strand']) < 0
                                             else "+"}
        return self.ids_proxy[db_id]

    def get_variant_id_from_protein_id(self,  **kwargs):
        """
        returns all information needed to instantiate a variation

        :param trans_id: A transcript ID (either ENSAMBLE (ENS) or RefSeq (NM, XN)
        :return: list of dicts -- containing all information needed for a variant initialization
        """
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        filter = None
        db_id = ""
        if "refseq" in kwargs:
            filter = "refseq_peptide"
            db_id = kwargs["refseq"]
        elif "ensemble" in kwargs:
            filter = "ensembl_peptide_id"
            db_id = kwargs["ensemble"]
        elif "swiss_accid" in kwargs:
            filter = "uniprot_swissprot_accession"
            db_id = kwargs["swiss_accid"]
        elif "swiss_gene" in kwargs:
            filter= "uniprot_swissprot"
            db_id = kwargs["swiss_gene"]
        else:
            warnings.war
        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(filter, str(db_id)) \
               + self.biomart_filter%("germ_line_variation_source", "dbSNP") \
               + self.biomart_attribute%("snp_ensembl_gene_id") \
               + self.biomart_attribute%("variation_name") \
               + self.biomart_attribute%("snp_chromosome_name") \
               + self.biomart_attribute%("chromosome_location") \
               + self.biomart_attribute%("allele") \
               + self.biomart_attribute%("snp_strand") \
               + self.biomart_attribute%("peptide_location") \
               + self.biomart_tail
        tsvreader = csv.DictReader(urllib2.urlopen(self.new_biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if not tsvselect:
            warnings.warn("No entry found for ID %s"%db_id)
            return None

        return tsvselect

    def get_variant_id_from_gene_id(self,  **kwargs):
        """
        returns all information needed to instantiate a variation

        :param trans_id: A transcript ID (either ENSAMBLE (ENS) or RefSeq (NM, XN)
        :return: list of dicts -- containing all information needed for a variant initialization
        """
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        filter = None
        db_id = ""
        if "hgnc_symbol" in kwargs:
            filter = "hgnc_symbol"
            db_id = kwargs["hgnc_symbol"]
        elif "ensemble" in kwargs:
            filter = "ensembl_gene_id"
            db_id = kwargs["ensemble"]
        elif "hgnc_id" in kwargs:
            filter = "hgnc_id"
            db_id = kwargs["hgnc_id"]
        elif "swiss_gene" in kwargs:
            filter = "uniprot_swissprot"
            db_id = kwargs["swiss_gene"]
        elif "external_gene_name" in kwargs:
            filter = "external_gene_name"
            db_id = kwargs["external_gene_name"]
        else:
            warnings.war
        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(filter, str(db_id)) \
               + self.biomart_filter%("germ_line_variation_source", "dbSNP") \
               + self.biomart_attribute%("snp_ensembl_gene_id") \
               + self.biomart_attribute%("variation_name") \
               + self.biomart_attribute%("snp_chromosome_name") \
               + self.biomart_attribute%("chromosome_location") \
               + self.biomart_attribute%("allele") \
               + self.biomart_attribute%("snp_strand") \
               + self.biomart_attribute%("peptide_location") \
               + self.biomart_tail
        tsvreader = csv.DictReader(urllib2.urlopen(self.new_biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        if not tsvselect:
            warnings.warn("No entry found for ID %s"%db_id)
            return None

        return tsvselect

    def get_protein_sequence_from_protein_id(self, **kwargs):
        """
        Returns the protein sequence for a given protein ID that can either be refeseq, uniprot or ensamble id

        :param kwargs:
        :return:
        """
        filter = None
        db_id = ""

        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        if "refseq" in kwargs:
            filter = "refseq_peptide"
            db_id = kwargs["refseq"]
        elif "ensemble" in kwargs:
            filter = "ensembl_peptide_id"
            db_id = kwargs["ensemble"]
        elif "swiss_accid" in kwargs:
            filter = "uniprot_swissprot_accession"
            db_id = kwargs["swiss_accid"]
        elif "swiss_gene" in kwargs:
            filter= "uniprot_swissprot"
            db_id = kwargs["swiss_gene"]
        else:
            warnings.warn("No correct transcript id")
            return None
        rq_n = self.biomart_head%(_db, _dataset) \
               + self.biomart_filter%(filter, str(db_id)) \
               + self.biomart_attribute%(filter) \
               + self.biomart_attribute%("peptide") \
               + self.biomart_tail
        print rq_n
        tsvreader = csv.DictReader(urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read().splitlines(), dialect='excel-tab')
        tsvselect = [x for x in tsvreader]
        print tsvselect
        if not tsvselect:
            warnings.warn("No entry found for ID %s"%db_id)
            return None

        return {EAdapterFields.SEQ: tsvselect[0]['Coding sequence'],
                                             EAdapterFields.GENE: tsvselect[0]['Associated Gene Name'],
                                             EAdapterFields.STRAND: "-" if int(tsvselect[0]['Strand']) < 0
                                             else "+"}

    #TODO: refactor ... function based on old code
    def get_all_variant_gene(self, locations,_db="hsapiens_gene_ensembl", _dataset='gene_ensembl_config'):
        """
        Fetches the important db ids and names for given chromosomal location
        :param chrom: integer value of the chromosome in question
        :param start: integer value of the variation start position on given chromosome
        :param stop: integer value of the variation stop position on given chromosome
        :return: The respective gene name, i.e. the first one reported
        """
        #TODO assert types
        #<!DOCTYPE Query><Query client="true" processor="TSV" limit="-1" header="1"><Dataset name="hsapiens_gene_ensembl" config="gene_ensembl_config"><Filter name="chromosomal_region" value="1:40367114:40367114,1:40702744:40702744,1:40705023:40705023,1:40771399:40771399,1:40777210:40777210,1:40881015:40881015,1:41235036:41235036,1:42048927:42048927,1:43002232:43002232,1:43308758:43308758,1:43393391:43630154,1:43772617:43772617,1:43772834:43772834" filter_list=""/><Attribute name="uniprot_genename"/></Dataset></Query>
        #queries = [self.biomart_filter%("uniprot_genename", ','.join(kwargs['genes'][x:x+300])) for x in xrange(0, len(kwargs['genes']), 300)]
        # if chrom + start + stop in self.gene_proxy:
        #     return self.gene_proxy[chrom + start + stop]
        #
        # rq_n = self.biomart_head \
        #     + self.biomart_filter%("chromosome_name", str(chrom))  \
        #     + self.biomart_filter%("start", str(start))  \
        #     + self.biomart_filter%("end", str(stop))  \
        #     + self.biomart_attribute%("uniprot_genename")  \
        #     + self.biomart_tail
        #
        # tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        # tsvselect = [x for x in tsvreader]
        # if tsvselect and tsvselect[0]:
        #     self.gene_proxy[chrom + start + stop] = tsvselect[0]['UniProt Gene Name']
        #     return tsvselect[0]['UniProt Gene Name']
        # else:
        #     logging.warn(','.join([str(chrom), str(start), str(stop)]) + ' does not denote a known gene location')
        #     return ''

    #TODO: refactor ... function based on old code
    def get_variant_ids(self, **kwargs):
        """
        Fetches the important db ids and names for given gene _or_ chromosomal location. The former is recommended.
        AResult is a list of dicts with either of the tree combinations:
            - 'Ensembl Gene ID', 'Ensembl Transcript ID', 'Ensembl Protein ID'
            - 'RefSeq Protein ID [e.g. NP_001005353]', 'RefSeq mRNA [e.g. NM_001195597]', first triplet
            - 'RefSeq Predicted Protein ID [e.g. XP_001720922]', 'RefSeq mRNA predicted [e.g. XM_001125684]', first triplet
        :keyword 'chrom': integer value of the chromosome in question
        :keyword 'start': integer value of the variation start position on given chromosome
        :keyword 'stop': integer value of the variation stop position on given chromosome
        :keyword 'gene': string value of the gene of variation
        :keyword 'transcript_id': string value of the gene of variation
        :return: The list of dicts of entries with transcript and protein ids (either NM+NP or XM+XP)
        """
        # TODO type assessment
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        ensemble_only = False
        query = None
        if len(kwargs) == 4 and 'chrom' in kwargs and 'start' in kwargs and 'stop' in kwargs and 'ensemble_only' in kwargs:
            ensemble_only = kwargs['ensemble_only']
            query = self.biomart_filter%("chromosome_name", kwargs['chrom'])  \
                + self.biomart_filter%("start", kwargs['start'])  \
                + self.biomart_filter%("end", kwargs['stop'])
        if len(kwargs) == 3 and 'chrom' in kwargs and 'start' in kwargs and 'stop' in kwargs:
            query = self.biomart_filter%("chromosome_name", kwargs['chrom'])  \
                + self.biomart_filter%("start", kwargs['start'])  \
                + self.biomart_filter%("end", kwargs['stop'])
        elif len(kwargs) == 2 and 'gene' in kwargs and 'ensemble_only' in kwargs:
            ensemble_only = kwargs['ensemble_only']
            query = self.biomart_filter%("uniprot_genename", kwargs['gene'])
        elif len(kwargs) == 2 and 'transcript_id' in kwargs and 'ensemble_only' in kwargs:
            ensemble_only = kwargs['ensemble_only']
            query = self.biomart_filter%("ensembl_transcript_id", kwargs['transcript_id'])
        elif len(kwargs) == 1 and 'gene' in kwargs:
            query = self.biomart_filter%("uniprot_genename", kwargs['gene'])
        else:
            warnings.warn("wrong arguments to get_variant_ids")

        rq_n = self.biomart_head%(_db, _dataset) \
            + query \
            + self.biomart_attribute%("uniprot_genename")  \
            + self.biomart_attribute%("ensembl_gene_id")  \
            + self.biomart_attribute%("ensembl_peptide_id")  \
            + self.biomart_attribute%("ensembl_transcript_id")  \
            + self.biomart_attribute%("strand")
        if not ensemble_only:
            rq_n += self.biomart_attribute%("refseq_mrna")  \
                + self.biomart_attribute%("refseq_peptide")
        rq_n += self.biomart_attribute%("uniprot_swissprot") + self.biomart_tail

        # logging.warning(rq_n)

        tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
        if ensemble_only:
            result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader}
        else:
            result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                  if x['RefSeq Protein ID [e.g. NP_001005353]'] and x['RefSeq mRNA [e.g. NM_001195597]']}

        if not ensemble_only:
            rq_x = self.biomart_head \
                + query  \
                + self.biomart_attribute%("uniprot_genename")  \
                + self.biomart_attribute%("ensembl_gene_id")  \
                + self.biomart_attribute%("ensembl_peptide_id")  \
                + self.biomart_attribute%("ensembl_transcript_id")  \
                + self.biomart_attribute%("refseq_peptide_predicted")  \
                + self.biomart_attribute%("refseq_mrna_predicted")  \
                + self.biomart_attribute%("strand")  \
                + self.biomart_tail

            tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_x)).read()).splitlines(), dialect='excel-tab')

            result2 = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                       if (x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]'])
                       or (not x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and not x['RefSeq mRNA predicted [e.g. XM_001125684]'])}
            result.update(result2)

        g = None
        for k, v in result.iteritems():
            if 'uniprot_swissprot' in v:
                g = v['uniprot_swissprot']
        self.ids_proxy[g] = result.values()
        return result.values()

    #TODO: refactor ... function based on old code
    def get_all_variant_ids(self,  **kwargs):
        """
        Fetches the important db ids and names for given gene _or_ chromosomal location. The former is recommended.
        AResult is a list of dicts with either of the tree combinations:
            - 'Ensembl Gene ID', 'Ensembl Transcript ID', 'Ensembl Protein ID'
            - 'RefSeq Protein ID [e.g. NP_001005353]', 'RefSeq mRNA [e.g. NM_001195597]', first triplet
            - 'RefSeq Predicted Protein ID [e.g. XP_001720922]', 'RefSeq mRNA predicted [e.g. XM_001125684]', first triplet
        :keyword 'locations': list of locations as triplets of integer values representing (chrom, start, stop)
        :keyword 'genes': list of genes as string value of the genes of variation
        :return: The list of dicts of entries with transcript and protein ids (either NM+NP or XM+XP)
        """
        _db = kwargs.get("_db","hsapiens_gene_ensembl")
        _dataset = kwargs.get("_dataset", "gene_ensembl_config")
        end_result = dict()
        # TODO type assessment
        ensemble_only = False
        query = None
        if 'ensemble_only' in kwargs:
            ensemble_only = kwargs['ensemble_only']
        if 'locations' in kwargs:
            pass
            #TODO
            # query = self.biomart_filter%("chromosome_name", kwargs['chrom'])  \
            #         + self.biomart_filter%("start", kwargs['start'])  \
            #         + self.biomart_filter%("end", kwargs['stop'])
        elif 'genes' in kwargs:
            queries = [self.biomart_filter%("uniprot_genename", ','.join(kwargs['genes'][x:x+250])) for x in xrange(0, len(kwargs['genes']), 250)]
            logging.warning('***'+self.biomart_filter%("uniprot_genename", ','.join(kwargs['genes'][x:x+250])) for x in xrange(0, len(kwargs['genes']), 250))
        else:
            logging.warning("wrong arguments to get_variant_ids")
        for query in queries:
            rq_n = self.biomart_head%(_db,_dataset) \
                + query \
                + self.biomart_attribute%("uniprot_genename")  \
                + self.biomart_attribute%("ensembl_gene_id")  \
                + self.biomart_attribute%("ensembl_peptide_id")  \
                + self.biomart_attribute%("ensembl_transcript_id")  \
                + self.biomart_attribute%("strand")
            if not ensemble_only:
                rq_n += self.biomart_attribute%("refseq_mrna")  \
                    + self.biomart_attribute%("refseq_peptide")
            rq_n += self.biomart_tail
            # rq_n += self.biomart_attribute%("uniprot_swissprot") + self.biomart_tail

            # logging.warning(rq_n)

            try:
                tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_n)).read()).splitlines(), dialect='excel-tab')
                if ensemble_only:
                    result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader}
                else:
                    result = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                          if x['RefSeq Protein ID [e.g. NP_001005353]'] and x['RefSeq mRNA [e.g. NM_001195597]']}
                end_result.update(result)
            except:
                logging.error('Bad Mart Query: '+rq_n)

            if not ensemble_only:
                rq_x = self.biomart_head%(_db, _dataset) \
                    + query  \
                    + self.biomart_attribute%("uniprot_genename")  \
                    + self.biomart_attribute%("ensembl_gene_id")  \
                    + self.biomart_attribute%("ensembl_peptide_id")  \
                    + self.biomart_attribute%("ensembl_transcript_id")  \
                    + self.biomart_attribute%("refseq_peptide_predicted")  \
                    + self.biomart_attribute%("refseq_mrna_predicted")  \
                    + self.biomart_attribute%("strand")  \
                    + self.biomart_tail

                try:
                    tsvreader = csv.DictReader((urllib2.urlopen(self.biomart_url+urllib2.quote(rq_x)).read()).splitlines(), dialect='excel-tab')

                    for x in tsvreader:
                        if (x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]']):
                            end_result.setdefault(x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID'], x)
                    # result2 = {x['Ensembl Gene ID']+x['Ensembl Transcript ID']+x['Ensembl Protein ID']: x for x in tsvreader
                    #            if (x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and x['RefSeq mRNA predicted [e.g. XM_001125684]'])
                    #     this line is a troublemaker and does not help       or (not x['RefSeq Predicted Protein ID [e.g. XP_001720922]'] and not x['RefSeq mRNA predicted [e.g. XM_001125684]'])}
                    # end_result.setdefault(result2)
                except:
                    logging.error('Bad Mart Query: '+rq_n)

        # g = None
        # for k, v in result.iteritems():
        #     if 'uniprot_swissprot' in v:
        #         g = v['uniprot_swissprot']
        # self.ids_proxy[g] = result.values()
        return end_result.values()
