from collections import defaultdict

__author__ = 'walzer'


from Fred2.Core.Generator import generate_peptides_from_variants
from Fred2.Core.Allele import Allele
from Fred2.Core.Peptide import Peptide
from Fred2.Core.Base import AEpitopePrediction
from Fred2.EpitopePrediction import EpitopePredictorFactory
from Fred2.IO.MartsAdapter import MartsAdapter
from Fred2.IO.EnsemblAdapter import EnsemblDB
import re
import csv
from Fred2.Core import Variant
from Fred2.Core import VariationType
from Fred2.Core import MutationSyntax
import xlrd
import re
rx = re.compile("[A-z]*[0-9]*")

filename_vars = "/home/walzer/bwSyncAndShare/LICA-FR/ssm_controlled.tsv"
filename_alleles = "/home/walzer/bwSyncAndShare/LICA-FR/2015-08-31 LiCa Typings_ML_final_MW_removedDateFormatting.xlsx"
filename_ensebl = "/home/walzer/dbs/Homo_sapiens.GRCh37.75.cds.all.fa"
xl_workbook = xlrd.open_workbook(filename_alleles)
xl_sheet = xl_workbook.sheet_by_index(0)

sect = None
typings = defaultdict(list)
for row_idx in range(1, xl_sheet.nrows):    # Iterate through rows
    dv = xl_sheet.row(row_idx)[1].value
    rs = rx.search(dv)
    donor = dv[rs.start():rs.end()]
    for i in range(2,6):
        hla_pref = 'HLA-A*' if i in [2,3] else 'HLA-B*'
        try:
            typings[donor].append(Allele(hla_pref+':'.join(str(xl_sheet.row(row_idx)[i].value).split(':')[0:2])))
        except Exception as e:
            if donor.startswith("CHC1079"):
                sect = xl_sheet.row(row_idx)
            continue


var_register = set()
all_vars = defaultdict(list)
all_peps9 = dict()

tmp_lst = list()
all_dict = list()


ma = MartsAdapter(biomart="http://grch37.ensembl.org")

with open(filename_vars, 'r') as f:
    dr = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE)
    for l in dr:
        #all_dict.append(l)
        #tmp_lst.append(l['consequence_type'])
        #rs = rx.search(l["icgc_sample_id"])
        rs = rx.search(l["submitted_sample_id"])
        #TODO for Variant: remember reference / build in a member for that! (in ICGC: assembly_version)
        isHomozygous = True if (l['tumour_genotype'].split('/')[0] == l['tumour_genotype'].split('/')[1]) else False
        if l['consequence_type'] == "missense_variant"\
                and l["icgc_mutation_id"] not in var_register:
            try:
                var = Variant(l["icgc_mutation_id"], VariationType.SNP, l['chromosome'], int(l['chromosome_start']),
                          l['reference_genome_allele'], l['mutated_to_allele'],
                          {l['transcript_affected']: MutationSyntax(l['transcript_affected'],
                                    int(filter(type(l['cds_mutation']).isdigit, l['cds_mutation'])),
                                    int(filter(type(l['aa_mutation']).isdigit, l['aa_mutation'])),
                                    "c."+l['cds_mutation'], "p."+l['aa_mutation'])},
                          isHomozygous, l['consequence_type'] == 'synonymous_variant')
                #all_vars[l["icgc_sample_id"][rs.start():rs.end()]].append(var)
                all_vars[l["submitted_sample_id"][rs.start():rs.end()]].append(var)
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
                    tr = ma.get_transcript_position(int(l['chromosome_start']), int(l['chromosome_start']),
                                                    l['gene_affected'], l['transcript_affected'])
                    t_pos = tr[0] if tr else '?'
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
                    tr = ma.get_transcript_position(int(l['chromosome_start']), int(l['chromosome_start']),
                                                    l['gene_affected'], l['transcript_affected'])
                    t_pos = str(tr[0]) if tr else '?'
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
# print ma.get_transcript_position('17953929', '17953943', 'ENSG00000074964', 'ENST00000361221')
# print ma.get_transcript_sequence('ENST00000361221')

# print len(all_vars), len(tmp_lst)
# from pygg import *
# import pandas
# df = pandas.DataFrame({'var_types': tmp_lst})
# #p = ggplot('data', aes('var_types')) + geom_histogram()
# #ggsave(p, "/tmp/out.png", data=df)
# from ggplot import *
# p = ggplot(aes('var_types'), data=df) + geom_histogram() + theme(axis_text_x=element_text(angle=90))
# f = p.draw()
# f.savefig('/tmp/t.png')


ep = EpitopePredictorFactory('netmhc')
#set !!! to be netMHC in $path dir hast to be first as well
# ln -s /home/walzer/immuno-tools/netMHC-3.4/netMHC-3.4 /home/walzer/immuno-tools/netMHC-3.4/netMHC

# test_allele = Allele("HLA-A*02:01")
# test_9ers = list()
# test_res = None
# test_gen_9ers = generate_peptides_from_variants(all_vars[all_vars.keys()[0]], 9, ma)
# c = 0
# for i in test_gen_9ers:
#     test_9ers.append(i)
#     c+=1
#     if c==3:
#         test_res = ep.predict(test_9ers, alleles=[test_allele])
#         break
#
# test_res
# test_res.iloc[0]
# peps = list(test_res.index)

ed = EnsemblDB()
ed.read_seqs(filename_ensebl)

for donor in all_vars:
    if not donor in typings:
        print donor, "not in allele list"
        continue
    all_9ers = generate_peptides_from_variants(all_vars[donor], 9, ed)
    ps = [p for p in all_9ers]
    donor_alleles = [a for a in typings[donor] if a.name in ep.supportedAlleles]
    try:
        all_peps9[donor] = ep.predict(all_9ers, alleles=donor_alleles)
    except Exception as e:
        print "Exception on", donor, "with so many peptides", len(ps)
        print e

for donor,df in all_peps9.iteritems():
    df.to_csv('/tmp/'+donor+'.csv')


