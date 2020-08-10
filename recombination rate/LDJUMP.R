#! /usr/bin/Rscript



#Before entring R
#module load R/3.6.0 R_packages/3.6.0
#

#library("devtools")
#devtools::install_github(repo="knausb/vcfR")
#install.packages("/home/venkat/bin/LDJump.tar.gz", repos=NULL, type="source")
library("LDJump")
library("Biostrings")
#VCF='batch.East.DPSCF300100.maf0.04.vcf.gz'
#ref_seq = "ref2.fa"
#chr = 100
#startofseq = 1
#endofseq = 559321


files <- read.table('/crex/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/flie_list.txt', header=F, col.names=c('Files'))
Autosomes <- read.table('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/Autosomes.txt', header=F, col.names=c('CHROM'))



for (val in files$Files) {
	tryCatch({
	Chrom=substr(val, 0, 11)
	out_file_name=paste('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDHat/out_put_rho/', val, '.table', sep="")
	if (file.exists(out_file_name)){}
	else {
	if (Chrom %in% Autosomes$CHROM) {
		setwd('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDHat/out_put_rho/')
		#LDJump(VCF, chr = chr, cores = 10, segLength = 10000, pathLDhat = '/home/venkat/bin/LDhat/', pathPhi = "/home/venkat/bin/PhiPack/", format = "vcf", refName = ref_seq, lengthofseq=559321, startofseq = startofseq, endofseq = endofseq)
		results <- LDJump(paste("/crex/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/windows/", val, sep=""), cores = 5, alpha = c(0.1, 0.05, 0.01), segLength = 10000, pathLDhat = '/home/venkat/bin/LDhat/', pathPhi = "/home/venkat/bin/PhiPack/Phi", format = "fasta",  start = 1, constant = F, status = T)
		final <- as.numeric(results[[2]])
		dnafile_1 <- read.fasta(file = paste("/crex/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDhelmet/windows/", val, sep=""))
		dnafile_1 <- dnafile_1[[1]]
		GC_final <- c(GC(dnafile_1[0:10000]), GC(dnafile_1[10000:20000]), GC(dnafile_1[20000:30000]), GC(dnafile_1[10000:40000]), GC(dnafile_1[10000:50000]))
		final_DB <- data.frame("GC" = GC_final, "r" = final)
		write.table(final_DB, paste('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Scripts/LDHat/out_put_rho/', val, '.table', sep=""))
		}
	}
	}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}

