import vcf
import ftplib
import os
import csv

# 1000 genomes data direct download via FTP in script
# ALI data available at https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000334.v1.p1

euros = []

panel = open("phase1_integrated_calls.20101123.ALL.panel.txt","r")

for row in panel:
    if row.split("\t")[2]=="EUR":
        euros += [row.split("\t")[0]]

panel.close()


# go through each chromosome
for i in range(1,23):
    ali_read=vcf.Reader(filename="ali_anno.vcf.gz")

    vcf1 = "ALL.chr" + str(i) + ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
    vcf2 = "ALL.chr" + str(i) + ".phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz.tbi"

    ftp = ftplib.FTP('ftp-trace.ncbi.nih.gov')
    ftp.login()
    ftp.cwd("1000genomes/ftp/release/20110521/")

    ftp.retrbinary("RETR " + vcf1, open(vcf1, 'wb').write)

    ftp = ftplib.FTP('ftp-trace.ncbi.nih.gov')
    ftp.login()
    ftp.cwd("1000genomes/ftp/release/20110521/")

    ftp.retrbinary("RETR " + vcf2, open(vcf2, 'wb').write)

    otg_read=vcf.Reader(filename=vcf1)

    SNPinfo = open("called_" + str(i) + "_snp_info.csv","wb")
    ALIdata = open("called_" + str(i) + "_ALI_data.csv","wb")
    OTGdata = open("called_" + str(i) + "_OTG_data.csv","wb")

    SNPwriter = csv.writer(SNPinfo,delimiter=",")
    SNPwriter.writerows([["SNP_POS","ALI_RS_Num","ALI_Sub_Num","OTG_RS_Num","OTG_Sub_Num","ALI_REF","ALI_ALT","OTG_REF","OTG_ALT"]])

    ALIwriter = csv.writer(ALIdata,delimiter=",")
    ALIwriter.writerows([["Subject","GT"]])

    OTGwriter = csv.writer(OTGdata,delimiter=",")
    OTGwriter.writerows([["Subject","GT"]])

    chrom_num = i

    cont1 = True

    while cont1:
        cur_ali = ali_read.next()
        if cur_ali.CHROM == str(chrom_num):
            cont1 = False

    cur_otg = otg_read.next()

    cont2 = True

    while cont2:

        if cur_ali.CHROM != str(chrom_num):
            cont2 = False

        if (cur_ali.CHROM==cur_otg.CHROM) and (cur_ali.POS==cur_otg.POS):

            cur_ali_samples = []
            cur_otg_samples = []

            pass_filter = True

            if ('QD' in cur_ali.INFO.keys()) and ('AB' in cur_ali.INFO.keys()) and ('SB' in cur_ali.INFO.keys()):
                if cur_ali.QUAL < 30.0:
                    pass_filter = False

                if cur_otg.QUAL < 30.0:
                    pass_filter = False

                if cur_ali.INFO['QD'] < 5.0:
                    pass_filter = False

                if cur_ali.INFO['AB'] > 0.75:
                    pass_filter = False

                if cur_ali.INFO['SB'] > -0.1:
                    pass_filter = False

            else:
                pass_filter = False




            total = 0
            nm = 0
            for sample in cur_ali.samples:
                total += 1
                if (sample.data.GT != './.') and (sample.data.GT != '.|.'):
                    nm += 1
                    cur_ali_samples += [[sample.sample]+[sample.data.GT]]

            if float(nm)/float(total) <= 0.9:
                pass_filter = False

            nm1 = 0
            for sample in cur_otg.samples:
                if ( (sample.data.GT != './.') and (sample.data.GT != '.|.') ) and (sample.sample in euros):
                    nm1 += 1
                    cur_otg_samples += [[sample.sample]+[sample.data.GT]]

            if float(nm1)/float(len(euros)) <= 0.9:
                pass_filter = False




            if pass_filter:
                ALIwriter.writerows(cur_ali_samples)
                OTGwriter.writerows(cur_otg_samples)

                SNPwriter.writerows([[cur_ali.POS,cur_ali.ID,len(cur_ali_samples),cur_otg.ID,len(cur_otg_samples),cur_ali.REF,cur_ali.ALT,cur_otg.REF,cur_otg.ALT]])

            cur_ali = ali_read.next()
            cur_otg = otg_read.next()

        elif (cur_ali.CHROM==cur_otg.CHROM) and (cur_ali.POS>cur_otg.POS):
            cur_otg = otg_read.next()
        elif (cur_ali.CHROM==cur_otg.CHROM) and (cur_ali.POS<cur_otg.POS):
            cur_ali = ali_read.next()

    SNPinfo.close()
    ALIdata.close()
    OTGdata.close()

    os.remove(vcf1)
    os.remove(vcf2)