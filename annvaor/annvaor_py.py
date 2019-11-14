# /usr/bin/env python
# coding=utf-8
###################################
#  Author : yunkeli
#  Version : 1.0(2016/3/1)
#  annovar unknown fileter
#  E-mail : 1316014512@qq.com
###################################
import json
import os,sys
import re
import thread
import threading
class Annovar2vcf():
	def __int__(self):
		global annotate_variation
		global annovar_database
		global convert2annovar
		global Gene_symbol_re
	def AnnovarPath(self,config):
		fileopen = open(config).read()
		return json.loads(fileopen)
	def AnnovarAll(self,vcffile,outdir,sample,annotate_variation,annovar_database,convert2annovar):
		c2a = convert2annovar+' -format vcf4 '+vcffile+'  > '+outdir+'/'+sample+'.avinput'
		print c2a
		c2g = convert2annovar+' -format vcf4 '+vcffile+'   -include -withzyg > '+outdir+'/'+sample+'.gffout'
		thousand_all = annotate_variation + ' -filter -dbtype 1000g2014oct_all -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput ' + annovar_database
		thousand_afr = annotate_variation + ' -filter -dbtype 1000g2014oct_afr -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput ' + annovar_database
		thousand_eas = annotate_variation + ' -filter -dbtype 1000g2014oct_eas -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput ' + annovar_database
		thousand_eur = annotate_variation + ' -filter -dbtype 1000g2014oct_eur -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput ' + annovar_database
		dbsnp = annotate_variation+' -filter -dbtype avsnp144  -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		ljb26_all  = annotate_variation+' -filter -dbtype ljb26_all  -otherinfo -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		ref  = annotate_variation+' -buildver hg19  -hgvs  -otherinfo -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		esp = annotate_variation+' -filter -dbtype esp6500siv2_all   -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		espaa = annotate_variation+' -filter -dbtype esp6500siv2_aa   -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		espea = annotate_variation+' -filter -dbtype esp6500siv2_ea   -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		cyb = annotate_variation+' -regionanno -dbtype cytoBand -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		cosmic_str = annotate_variation+' -filter -dbtype cosmic70   -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		clinvar_str = annotate_variation+' -filter -dbtype clinvar_20160302   -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		nci60_str = annotate_variation+' -filter -dbtype nci60   -buildver hg19 -out '+outdir+'/'+sample+' ' + outdir+'/'+sample+'.avinput '+annovar_database
		phastConsElements46way_str = annotate_variation+'  -regionanno -build hg19 -out '+outdir+'/'+sample+' -dbtype phastConsElements46way ' + outdir+'/'+sample+'.avinput '+annovar_database
		tfbsConsSites_str = annotate_variation+'  -regionanno -build hg19 -out '+outdir+'/'+sample+' -dbtype tfbsConsSites ' + outdir+'/'+sample+'.avinput '+annovar_database
		wgRna_str = annotate_variation+' -regionanno -build hg19 -out '+outdir+'/'+sample+' -dbtype wgRna ' + outdir+'/'+sample+'.avinput '+annovar_database
		targetScanS_str = annotate_variation+'  -regionanno -build hg19 -out '+outdir+'/'+sample+' -dbtype targetScanS ' + outdir+'/'+sample+'.avinput '+annovar_database
		gwasCatalog_str = annotate_variation+'  -regionanno -build hg19 -out '+outdir+'/'+sample+' -dbtype gwasCatalog ' + outdir+'/'+sample+'.avinput '+annovar_database
		
		# print c2g
		# print thousand_afr
		# print wgRna_str
		# sys.exit()
		c2ainfo = os.system(c2a)
		alllist = []
		def c2gtd():
			c2ginfo = os.system(c2g)
			alllist.append("c2g")
		def cosmictd():
			cosmicinfo = os.system(cosmic_str)
			alllist.append("cosmic70")
		def clinvartd():
			clinvarfo = os.system(clinvar_str)
			alllist.append("clinvar")
		def nci60td():
			nciinfo = os.system(nci60_str)
			alllist.append("nci60")
		def t_afr():
			t_allinfo = os.system(thousand_afr)
			alllist.append("thousand_afr")
		def t_eas():
			t_allinfo = os.system(thousand_eas)
			alllist.append("eas")
		def t_eur():
			t_allinfo = os.system(thousand_eur)
			alllist.append("thousand_eur")
		def t_alltd():
			t_allinfo = os.system(thousand_all)
			alllist.append("thousand_all")
		def dstd():
			dsinfo = os.system(dbsnp)
			alllist.append("dbsnp")
		def ljb26_alltd():
			ljb26_allinfo = os.system(ljb26_all)
			alllist.append("ljb26_all")
		def reftd():
			refinfo = os.system(ref)
			alllist.append("ref")
		def esptd():
			espinfo = os.system(esp)
			alllist.append('esp6500')
		def espaatd():
			espinfo = os.system(espaa)
			alllist.append('esp6500aa')
		def espeatd():
			espinfo = os.system(espea)
			alllist.append('esp6500ea')
		def cybtd():
			cybinfo = os.system(cyb)
			alllist.append('cytoBand')
		def phast():
			cybinfo = os.system(phastConsElements46way_str)
			alllist.append('phastConsElements46way_str')
		def tfbs():
			cybinfo = os.system(tfbsConsSites_str)
			alllist.append('tfbsConsSites_str')
		def wgRna():
			cybinfo = os.system(wgRna_str)
			alllist.append(wgRna)
		def target():
			cybinfo = os.system(targetScanS_str)
			alllist.append('target')
		def gwas():
			cybinfo = os.system(gwasCatalog_str)
			alllist.append('gwas')
		if c2ainfo == 0:
			t1 = threading.Thread(target=c2gtd)
			t2 = threading.Thread(target=t_alltd)
			t3 = threading.Thread(target=dstd)
			t4 = threading.Thread(target=ljb26_alltd)
			t5 = threading.Thread(target=reftd)
			t6 = threading.Thread(target=esptd)
			t7 = threading.Thread(target=cybtd)
			t8 = threading.Thread(target=cosmictd)
			t9 = threading.Thread(target=clinvartd)
			t10 = threading.Thread(target=nci60td)
			t11 = threading.Thread(target=t_eas)
			t12 = threading.Thread(target=t_eur)
			t15 = threading.Thread(target=t_afr)
			t13 = threading.Thread(target=espaatd)
			t14 = threading.Thread(target=espeatd)
			t16 = threading.Thread(target=phast)
			t17 = threading.Thread(target=tfbs)
			t18 = threading.Thread(target=wgRna)
			t19 = threading.Thread(target=target)
			t20 = threading.Thread(target=gwas)
			t1.start()
			t2.start()
			t3.start()
			t4.start()
			t5.start()
			t6.start()
			t7.start()
			t8.start()
			t9.start()
			t10.start()
			t11.start()
			t12.start()
			t13.start()
			t14.start()
			t15.start()
			t16.start()
			t17.start()
			t18.start()
			t19.start()
			t20.start()
			while True:
				if len(alllist) == 20:
					#print alllist
					return 0
		else:
			print 'annovar is erro please check log'
			sys.exit()
	def Gene_symbol_re(self,content):
		link = re.compile("\(.*?\)")
		info = re.sub(link,"",content)
		return info
	# def inhouse(self,data):
		# fileopen = open(data).readlines()
		# inhousedict = {}
		# for line in fileopen:
			# linelist = line.split()
			# key = linelist[0]+'|'+linelist[1]+'|'+linelist[2]+'|'+linelist[4]
			# value = linelist[-1]
			# print key
			# inhousedict[key] = value
		# return inhousedict
	def dictavsnp144(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_avsnp144_dropped').readlines()
		snp144dict = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			snp144dict[key] = values
		# fileopen.close()
		return snp144dict
	def dict1000gall(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_ALL.sites.2014_10_dropped').readlines()
		all1000gdict = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			all1000gdict[key] = values
		# fileopen.close()
		return all1000gdict
	def dict1000gafr(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_AFR.sites.2014_10_dropped').readlines()
		afr1000gdict = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			afr1000gdict[key] = values
		# fileopen.close()
		return afr1000gdict
	def dict1000geur(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_EUR.sites.2014_10_dropped').readlines()
		eur1000gdict = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			eur1000gdict[key] = values
		# fileopen.close()
		return eur1000gdict
	def dict1000geas(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_EAS.sites.2014_10_dropped').readlines()
		eas1000gdict = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			eas1000gdict[key] = values
		# fileopen.close()
		return eas1000gdict
	def dictesp(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_esp6500siv2_all_dropped').readlines()
		espdict = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			espdict[key] = values
		# fileopen.close()
		return espdict
	def dictespaa(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_esp6500siv2_aa_dropped').readlines()
		espdictaa = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			espdictaa[key] = values
		# fileopen.close()
		return espdictaa
	def dictespea(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_esp6500siv2_ea_dropped').readlines()
		espdictea = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			espdictea[key] = values
		# fileopen.close()
		return espdictea
	def dictcyb(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_cytoBand').readlines()
		exac03dict = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			exac03dict[key] = values
		return exac03dict
	def dictvariant(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.variant_function').readlines()
		variantdict = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = [linelist[0],linelist[1]]
			variantdict[key] = values
		# fileopen.close()
		return variantdict
	def dictexonic(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.exonic_variant_function').readlines()
		exonicdict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[3]+'|'+linelist[4]+'|'+linelist[5]+'|'+linelist[7]
			values = [linelist[1].split()[0],linelist[2]]
			exonicdict[key] = values
		# fileopen.close()
		return exonicdict
	def dictljb(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_ljb26_all_dropped').readlines()
		ljbdict = {}
		for line in fileopen:
			linelist = line.split()
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			v = linelist[1].split(',')
			if len(v) == 25:
				ljbdict[key] = v
		return ljbdict
	def dictclinvar(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_clinvar_20160302_dropped').readlines()
		clinvardict = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			clinvardict[key] = values
		return clinvardict
	def cosmic70(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_cosmic70_dropped').readlines()
		cosmic = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			cosmic[key] = values
		return cosmic
	def nci60(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_nci60_dropped').readlines()
		nci = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			nci[key] = values
		return nci
	def phastdict(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_phastConsElements46way').readlines()
		phast = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			phast[key] = values
		return phast
	def tfbsdict(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_tfbsConsSites').readlines()
		tfbs = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			tfbs[key] = values
		return tfbs
	def wgRnadict(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_wgRna').readlines()
		wgRna = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			wgRna[key] = values
		return wgRna
	def targetdict(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_targetScanS').readlines()
		target = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			target[key] = values
		return target
	def gwasdict(self,outdir,sample):
		fileopen = open(outdir+'/'+sample+'.hg19_gwasCatalog').readlines()
		gwas = {}
		for line in fileopen:
			linelist = line.split('\t')
			key = linelist[2]+'|'+linelist[3]+'|'+linelist[4]+'|'+linelist[6]
			values = linelist[1]
			gwas[key] = values
		return gwas
	def MergeAll(self,outdir,sample,exonicdict,snp144dict,all1000gdict,cybdict,espdict,variantdict,ljbdict,afr1000gdict,eas1000gdict,eur1000gdict,espdictaa,espdictea,clinvar,cosmic,nci,phast,tfbs,wgRna,target,gwas):
		def Gene_symbol_re(content):
			link = re.compile("\(.*?\)")
			info = re.sub(link,"",content)
			return info
		fileopen9 = open(outdir+'/'+sample+'.gffout').readlines()
		head = '#Chr\tStart\tEnd\tRef\tAlt\tGenotype\tref/alt(dpth)\tFilter\tquality_score\tdepth\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tphastConsElements46way\ttfbsConsSites\tcytoBand\twgRna\ttargetScanS\tgwasCatalog\t1000g2014oct_all\t1000g2014oct_afr\t1000g2014oct_eas\t1000g2014oct_eur\tavsnp144\tGenotype\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\tPolyphen2_HVAR_score\tPolyphen2_HVAR_pred\tLRT_score\tLRT_pred\tMutationTaster_score\tMutationTaster_pred\tMutationAssessor_score\tMutationAssessor_pred\tFATHMM_score\tFATHMM_pred\tRadialSVM_score\tRadialSVM_pred\tLR_score\tLR_pred\tVEST3_score\tCADD_raw\tCADD_phred\tGERP++_RS\tphyloP46way_placental\tphyloP100way_vertebrate\tSiPhy_29way_logOdds\tesp6500si_all\tesp6500si_ea\tesp6500si_aa\tclinvar_20160802\tcosmic70\tnci60\n'
		filesave = open(outdir+'/'+sample+'.gff','w+')
		#filesave_stat = open(outdir+'/'+sample+'.stat','w+')
		filesave.write(head)
		#print head
		
		count = 0
		unknown_count = 0
		for x in fileopen9:
			p =  x.split()
			key = p[0]+'|'+p[1]+'|'+p[2]+'|'+p[4]
			pend = p[0:6]
			#pend_unknown = p[0:8]
			qc = p[-1].split(':')
			if len(qc) == 5:
				pend.append(qc[1])
				#pend_unknown.append(qc[1])
			else:
				pend.append('-')
				#pend_unknown.append('-')
			pend.append(p[14])
			pend.extend(p[6:8])
			if key in variantdict:
				pend.append(variantdict[key][0])
			else:
				pend.append('-')
			if key in variantdict:
				pend.append(Gene_symbol_re(variantdict[key][1]))
			else:
				pend.append('-')
			if key in variantdict:
				pend.append(variantdict[key][1])
			else:
				pend.append('-')
				
			if key in exonicdict:
			#	print key
				pend.append(exonicdict[key][0])
				pend.append(exonicdict[key][1])
			else:
				pend.append('unknown')
				pend.append('-')
			if key in phast:
				pend.append(phast[key])
			else:
				pend.append('-')
			if key in tfbs:
				pend.append(tfbs[key])
			else:
				pend.append('-')				
			if key in cybdict:
				pend.append(cybdict[key])
			else:
				pend.append('-')
			if key in wgRna:
				pend.append(wgRna[key])
			else:
				pend.append('-')
			if key in target:
				pend.append(target[key])
			else:
				pend.append('-')
			if key in gwas:
				pend.append(gwas[key])
			else:
				pend.append('-')
			if key in all1000gdict:
				pend.append(all1000gdict[key])
			else:
				pend.append('-')
			if key in afr1000gdict:
				pend.append(afr1000gdict[key])
			else:
				pend.append('-')
					
			if key in eas1000gdict:
				pend.append(eas1000gdict[key])
			else:
				pend.append('-')
				
			if key in eur1000gdict:
				pend.append(eur1000gdict[key])
			else:
				pend.append('-')				
			if key in snp144dict:
				pend.append(snp144dict[key])
			else:
				pend.append('-')
			if '.' not in p[3] and '.' not in p[4] and '-' not in p[3:5] and len(p[3]) == 1 and len(p[4]) == 1 :
				if p[5].strip() == "hom":
					pend.append(p[4]+p[4])
				else:
					pend.append(p[3]+p[4])
			else:
				pend.append('-')
				# if key in 
			if key in ljbdict:
				pend = pend+ljbdict[key]
				# print pend
			else:
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
				pend.append('-')
			if key in espdict:
				pend.append(espdict[key])
			else:
				pend.append('-')
			if key in espdictea:
				pend.append(espdictea[key])
			else:
				pend.append('-')
			if key in espdictaa:
				pend.append(espdictaa[key])
			else:
				pend.append('-')
			
			if key in clinvar:
				pend.append(clinvar[key])
			else:
				pend.append('-')
			if key in cosmic:
				pend.append(cosmic[key])
			else:
				pend.append('-')
			if key in nci:
				pend.append(nci[key])
			else:
				pend.append('-')
			# print pend
			
			filesave.write("\t".join(pend)+"\n")
		filesave.close()
if  __name__ == "__main__":
	import argparse
	usage='annovar_py.py [--outdir] [--inputvcf] [--samplename] [--config]'
	if len(sys.argv) < 2:
		print usage
		print 'please cmd "python annovar_py.py -h "'
		sys.exit()
	else:
		p = argparse.ArgumentParser(usage='annovar_py.py [--outdir] [--inputvcf] [--samplename] [--config]', description="annovar vcf to gff")
		p.add_argument('-i','--inputvcf',   type=str, help = 'input   dir eg: ./sample.vcf',required=True)
		p.add_argument('-s','--samplename', type=str, help = 'sample name',required=True)
		p.add_argument('-o','--outdir',     type=str, help ='oss document directory eg: ./yunkeli',required=True)
		p.add_argument('-c','--config', type=str, help ='config file',required=True)
		args = p.parse_args()
		inputvcf = args.inputvcf
		samplename = args.samplename
		outdir = args.outdir
		config = args.config
		project = Annovar2vcf()
		path = project.AnnovarPath(config)
		print path
		annotate_variation = path["annotate_variation"]
		annovar_database = path["annovar_database"]
		convert2annovar = path["convert2annovar"]
		otherdatabase = path["otherdatabase"]
		stat = project.AnnovarAll(inputvcf,outdir,samplename,annotate_variation,annovar_database,convert2annovar)
		# print stat
		snpdict = project.dictavsnp144(outdir,samplename)
		all1000dict = project.dict1000gall(outdir,samplename)
		ljbdict = project.dictljb(outdir,samplename)
		cybdict = project.dictcyb(outdir,samplename)
		variantdict	= project.dictvariant(outdir,samplename)
		dict1000gafr	= project.dict1000gafr(outdir,samplename)
		dict1000geas	= project.dict1000geas(outdir,samplename)
		dict1000geur	= project.dict1000geur(outdir,samplename)
		dictespaa	= project.dictespaa(outdir,samplename)
		dictespea	= project.dictespea(outdir,samplename)
		dict1clinvar	= project.dictclinvar(outdir,samplename)
		cosmic70	= project.cosmic70(outdir,samplename)
		nci60	= project.nci60(outdir,samplename)
		phast	= project.phastdict(outdir,samplename)
		tfbs	= project.tfbsdict(outdir,samplename)
		wgRna	= project.wgRnadict(outdir,samplename)
		target	= project.targetdict(outdir,samplename)
		gwas	= project.gwasdict(outdir,samplename)
		#print variantdict
		exonicdict = project.dictexonic(outdir,samplename)
		espdict = project.dictesp(outdir,samplename)
		annotate_m = project.MergeAll(outdir,samplename,exonicdict,snpdict,all1000dict,cybdict,espdict,variantdict,ljbdict,dict1000gafr,dict1000geas,dict1000geur,dictespaa,dictespea,dict1clinvar,cosmic70,nci60,phast,tfbs,wgRna,target,gwas)
		#print inhousedict
