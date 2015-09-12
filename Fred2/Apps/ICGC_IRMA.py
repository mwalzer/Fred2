from collections import defaultdict

__author__ = 'walzer'


from Fred2.IO import FileReader
from Fred2.Core.Generator import generate_peptides_from_protein
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import AEpitopePrediction
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.IO.MartsAdapter import MartsAdapter
import re
import csv
from Fred2.Core import Variant
from Fred2.Core import VariationType
from Fred2.Core import MutationSyntax

filename = "/home/walzer/bwSyncAndShare/LICA-FR/ssm_controlled.tsv"
rx = re.compile("[A-z]*[0-9]*")
isHomozygous = True #TODO find out via 'tumour_genotype'
var_register = set()
all_vars = defaultdict(list)

tmp_lst = list()
all_dict = list()

with open(filename, 'r') as f:
    dr = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    for l in dr:
        all_dict.append(l)
        tmp_lst.append(l['consequence_type'])
        rs = rx.search(l["icgc_sample_id"])
        #TODO for Variant: remember reference / build in a member for that! (in ICGC: assembly_version)
        if l['consequence_type'] == "missense_variant"\
                and l["icgc_mutation_id"] not in var_register:
            try:
                var = Variant(l["icgc_mutation_id"], VariationType.SNP, l['chromosome'], int(l['chromosome_start']),
                          l['reference_genome_allele'], l['mutated_to_allele'],
                          {l['transcript_affected']: MutationSyntax(l['transcript_affected'],
                                    int(filter(type(l['cds_mutation']).isdigit, l['cds_mutation'])),
                                    int(filter(type(l['aa_mutation']).isdigit, l['aa_mutation'])),
                                    "c."+l['cds_mutation'], "p."+l['aa_mutation'])}, isHomozygous, l['consequence_type'] == 'synonymous_variant')
                all_vars[l["icgc_sample_id"][rs.start():rs.end()]].append(var)
                var_register.add(l["icgc_mutation_id"])
            except Exception as e:
                print "ms"
                print l
                print e
                break
        elif l['consequence_type'] == "frameshift_variant"\
                and l["icgc_mutation_id"] not in var_register:
            try:
                if l["mutation_type"].startswith("insertion"):
                    vt = VariationType.FSINS
                else:
                    vt = VariationType.FSDEL
                if l['cds_mutation']:
                    t_pos = int(filter(type(l['cds_mutation']).isdigit, l['cds_mutation']))
                    t_syn = l['cds_mutation']
                    if not t_syn.startswith("c."):
                        t_syn = "c." + t_syn
                else:
                    t_pos = '?'
                    t_syn = "c.?"
                var = Variant(l["icgc_mutation_id"], vt, l['chromosome'], int(l['chromosome_start']),
                          l['reference_genome_allele'].strip('-'), l['mutated_to_allele'].strip('-'),
                          {l['transcript_affected']: MutationSyntax(l['transcript_affected'],
                                    t_pos,
                                    int(filter(type(l['aa_mutation']).isdigit, l['aa_mutation'])),
                                    t_syn,
                                    "p."+l['aa_mutation'])},
                          isHomozygous, False)
                all_vars[l["icgc_sample_id"][rs.start():rs.end()]].append(var)
                var_register.add(l["icgc_mutation_id"])
            except Exception as e:
                print "fs"
                print l
                print e
                break
        elif l['consequence_type'].startswith("inframe_")\
                and l["icgc_mutation_id"] not in var_register:
            try:
                if l["mutation_type"].endswith("insertion"):
                    vt = VariationType.INS
                else:
                    vt = VariationType.DEL
                if l['cds_mutation']:
                    t_pos = int(filter(type(l['cds_mutation']).isdigit, l['cds_mutation']))
                    t_syn = l['cds_mutation']
                    if not t_syn.startswith("c."):
                        t_syn = "c." + t_syn
                else:
                    t_pos = '?'
                    t_syn = "c.?"
                var = Variant(l["icgc_mutation_id"], vt, l['chromosome'], int(l['chromosome_start']),
                          l['reference_genome_allele'].strip('-'), l['mutated_to_allele'].strip('-'),
                          {l['transcript_affected']: MutationSyntax(l['transcript_affected'],
                                    t_pos,
                                    int(filter(type(l['aa_mutation']).isdigit, l['aa_mutation'])),
                                    t_syn,
                                    "p."+l['aa_mutation'])},
                          isHomozygous, False)
                all_vars[l["icgc_sample_id"][rs.start():rs.end()]].append(var)
                var_register.add(l["icgc_mutation_id"])
            except Exception as e:
                print "indel"
                print l
                print e
                break

#TODO what about exon_loss_variant and intragenic_variant s?

print len(all_vars), len(tmp_lst)
from pygg import *
import pandas
df = pandas.DataFrame({'var_types':tmp_lst})
#p = ggplot('data', aes('var_types')) + geom_histogram()
#ggsave(p, "/tmp/out.png", data=df)

from ggplot import *
p = ggplot(aes('var_types'), data=df) + geom_histogram() + theme(axis_text_x=element_text(angle=45))
f = p.draw()
f.savefig('/tmp/t.png')

#TODO read ensembl

#TODO problem in ICGC: no cds_mutation entry about the transcript position

#TODO build transcripts
#TODO build proteins
#TODO make predictions