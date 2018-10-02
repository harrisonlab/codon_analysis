# Extract gene IDs in list form belongs to High/Low expressed genes across conditions
# This function requies a table consisting of 10 columns
# Column1: gene ID, column2: gene Length, Columns 3 to 10: Average FPKM values per condition per column
extractID <- function(mat) {
	high <- list()
	low <- list()
	ids <- list()
	for(i in 3:10) {
		high[[i]] <- subset(matH <- mat[,c(1,2,i)],matH[,3] >= 148)
	}
	for(i in 3:10) {
		low[[i]] <- subset(matL <- mat[,c(1,2,i)], matL[,3] >= 5 & matL[,3] <= 10)
	}
	names(high) <- c("NA","NA","RH1","RH2","RH3","RH4","RH5","RH6","RH7","RH8")
	names(low) <- c("NA","NA","RH1","RH2","RH3","RH4","RH5","RH6","RH7","RH8")
	ids[[1]] <- high
	ids[[2]] <- low
	names(ids) <- c("High","Low")
	return(ids)
}
# End of extractID #
#
#
#
# Extract gene sequences for list of gene IDs for RSCU task
# This function requires a fasta formatted sequence file and IDs extracted using "extractID" function
extractSeq <- function(ids,fasta) {
	Seq <- list()
	highSeq <- list()
	lowSeq <- list()
	for(i in 3:10) {
		highSeq[[i]] <- fasta[ids[[1]][[i]][[1]]]
	}
	for(i in 3:10) {
		lowSeq[[i]] <- fasta[ids[[2]][[i]][[1]]]
	}
	names(highSeq) <- c("NA","NA","RH1","RH2","RH3","RH4","RH5","RH6","RH7","RH8")
	names(lowSeq) <- c("NA","NA","RH1","RH2","RH3","RH4","RH5","RH6","RH7","RH8")
	Seq[[1]] <- highSeq
	Seq[[2]] <- lowSeq
	names(Seq) <- c("High","Low")
	return(Seq)
}
# End of extractSeq #
#
#
#
# calRSCU function calculates RSCU values for each of the codon of each sequence and creates a matrix
# Usage: testMat <- calRSCU("test.fa")
#calRSCU <- function(fasta) {
#	#Seq <- readDNAStringSet(fasta)
#	Seq <- fasta
#	SeqMat <- c()
#	for(i in 1:length(Seq)) {
#		SeqMat <- cbind(SeqMat, SADEG.RSCU(Seq[i]))
#	}
#	SeqMat <- t(SeqMat)
#	rownames(SeqMat) <- names(Seq)
#	SeqMat <- SeqMat[,-c(38,55,62,63,64)]
#	return(SeqMat)
#}
# End of calRSCU #
#
#
###
###
###
processMat <- function(matH,matL,codon,codonw) {
	fus <- list()
	for(i in 3:10) {
		matH1 <- matrix(round(unlist(matH[[i]]),3),ncol=64,byrow=T)
		colnames(matH1) <- names(matH[[i]][[1]])
		rownames(matH1) <- names(matH[[i]])
		print(dim(matH1))
		print("------------------------------------")
		matH1 <- matH1[,-c(38,55,62,63,64)]

		matL1 <- matrix(round(unlist(matL[[i]]),3),ncol=64,byrow=T)
		colnames(matL1) <- names(matL[[i]][[1]])
		rownames(matL1) <- names(matL[[i]])
		print(dim(matL1))
		print("------------------------------------")
		matL1 <- matL1[,-c(38,55,62,63,64)]

		print(dim(matH1))
		print(dim(matL1))

		#print(head(matH1))
		#print("------------------------------------")
		#print(head(matL1))
		#print("------------------------------------")
		#print(codon)
		#print(dim(codonw))

		fus[[i]] <- deltaRSCU(testH=matH1,testL=matL1,codonw=codonw,codon=codon)


	}
	names(fus) <- c("NA","NA","RH1","RH2","RH3","RH4","RH5","RH6","RH7","RH8")
	print(fus)

	capture.output(fus, file = "fus.txt")
}
# End of processMat #
###
###
###
# deltaRSCU
# This function performs t-test between RSCU values of highly and lowly expressed genes of a condition
deltaRSCU <- function (testH,testL,codonw,codon) {

	tmp <- list()
	for(i in 1:dim(testL)[2]) {
		tmpL <- as.numeric(testL[1:dim(testL)[1],i])
		tmpH <- as.numeric(testH[1:dim(testH)[1],i])
		tmp[[i]] <- t.test(tmpH,tmpL,var.equal = T)
	}
	print(tmp[[1]])

	#make matrix 3 wide and 59 long put mean of x mean of y and p.value along with the rowname as a codon
	fusTest <- matrix(data = NA, nrow = 59, ncol = 4, byrow = FALSE,dimnames = NULL)

	for (i in 1:59) {
		fusTest[i,1]= as.numeric(tmp[[i]]$estimate[1])
		fusTest[i,2]=as.numeric(tmp[[i]]$estimate[2])
		fusTest[i,3]=tmp[[i]]$p.value
	}
	print(head(fusTest))

	rownames(fusTest)=codon

	print(codon)

	print(dim(fusTest))

	for (i in 1:59) {
		if (fusTest[i,1]>fusTest[i,2] & fusTest[i,3]<1e-2){
			fusTest[i,4]=3
		}
		if (fusTest[i,1]<fusTest[i,2] & fusTest[i,3]<1e-2){
			fusTest[i,4]=1
		}
		if (fusTest[i,1]<fusTest[i,2] & fusTest[i,3]>1e-2){
			fusTest[i,4]=2
		}
		if (fusTest[i,1]>fusTest[i,2] & fusTest[i,3]>1e-2){
			fusTest[i,4]=2
		}
	}
	print(head(fusTest))

	arrangedfus = matrix(data = NA, nrow = 64, ncol = 2, byrow = FALSE, dimnames = NULL)

	for (i in 1:64) {
		if ( length(which (rownames(fusTest)==codonw[i,]) )<1  ) {
			arrangedfus[i,1]=codonw[i,]
			arrangedfus[i,2]=2
		}
		else {
			arrangedfus[i,1]=rownames(fusTest)[which (rownames(fusTest)==codonw[i,])]
			arrangedfus[i,2]=as.numeric(fusTest[ which (rownames(fusTest)==codonw[i,]),4])
		}
	}
	print(arrangedfus)
	fusCOA =matrix(as.numeric(arrangedfus[,2]), nrow = 4, ncol = 16, byrow = TRUE, dimnames = NULL)
	print("Writing results to files!")

	#write.table(data.frame(arrangedfus), file = "fus_arr_new.txt",append = FALSE, sep = ",",row.names = FALSE, col.names = FALSE,qmethod = "escape")
	#write.table(fusCOA, file = "fus_new.coa",append = FALSE, sep = ",",row.names = FALSE, col.names = FALSE)

	return(fusCOA)
}
# End of deltaRSCU #
#
#
# No need of calculateRSCU function. Directly run the statements twice!
#calculateRSCU <- function(testSeq) {
#	matH <- lapply(testSeq[[1]], function(x) {lapply(x,SADEG.RSCU)})
#	matL <- lapply(testSeq[[2]], function(x) {lapply(x,SADEG.RSCU)})
#}
# End of calculateRSCU #
###
###
###

processMat2 <- function(matH,matL,codonw,Seq) {
        fus <- list()
        for(i in 3:10) {
                matH1 <- matrix(round(unlist(matH[[i]]),3),ncol=64,byrow=T)
                colnames(matH1) <- names(matH[[i]][[1]])
                rownames(matH1) <- names(matH[[i]])
                print(dim(matH1))
                print("------------------------------------")
                matH1 <- matH1[,-c(38,55,62,63,64)]

                matL1 <- matrix(round(unlist(matL[[i]]),3),ncol=64,byrow=T)
                colnames(matL1) <- names(matL[[i]][[1]])
                rownames(matL1) <- names(matL[[i]])
                print(dim(matL1))
                print("------------------------------------")
                matL1 <- matL1[,-c(38,55,62,63,64)]

                print(dim(matH1))
                print(dim(matL1))

                print(head(matH1))
                print("------------------------------------")
                print(head(matL1))
                print("------------------------------------")
		codon <- colnames(matL1)
                print(codon)
                print(dim(codonw))

                fus[[i]] <- deltaRSCU2(testH=matH1,testL=matL1,codonw=codonw,codon=codon)


        }
        #names(fus) <- c("NA","NA","RH1","RH2","RH3","RH4","RH5","RH6","RH7","RH8")
        #print(fus)

        #capture.output(as.matrix(fus[[3]]), file = "fus.txt")

	for(i in 4:4) {
		print("Running codonW")
		write.table(fus[[i]], file="fus.coa",append = FALSE, sep = ",",row.names = FALSE, col.names = FALSE)
		fasta <- c(Seq[[1]][[i]],Seq[[2]][[i]])
		writeXStringSet(fasta,"input.dat")

		system("/home/sobczm/popgen/codon/codonW/codonw input.dat -all_indices -nomenu -silent -fop_file fus.coa")
	}

	#write.table(fus[[3]], file = "fus.coa",append = FALSE, sep = ",",row.names = FALSE, col.names = FALSE)

	#writeXStringSet(Seq[[1]][[4]],"seq.fasta")

	#print(fus[3])
	#print("Running codonW")
	#system("/home/sobczm/popgen/codon/codonW/codonw seq.fasta -all_indices -nomenu -silent -fop_file fus.coa")
}

deltaRSCU2 <- function (testH,testL,codonw,codon) {

        tmp <- list()
        for(i in 1:dim(testL)[2]) {
                tmpL <- as.numeric(testL[1:dim(testL)[1],i])
                tmpH <- as.numeric(testH[1:dim(testH)[1],i])
                tmp[[i]] <- t.test(tmpH,tmpL,var.equal = T)
        }
        print(tmp[[1]])

        #make matrix 3 wide and 59 long put mean of x mean of y and p.value along with the rowname as a codon
        fusTest <- matrix(data = NA, nrow = 59, ncol = 4, byrow = FALSE,dimnames = NULL)

        for (i in 1:59) {
                fusTest[i,1]= as.numeric(tmp[[i]]$estimate[1])
                fusTest[i,2]=as.numeric(tmp[[i]]$estimate[2])
                fusTest[i,3]=tmp[[i]]$p.value
        }
	print("*****************************")
        print(head(fusTest))
	print("*****************************")
	print(codon)
	print(length(codon))
	print(rownames(fusTest))
	print(dim(fusTest))
        rownames(fusTest)=codon

	print(codon)
	print("*****************************")

        for (i in 1:59) {
                if (fusTest[i,1]>fusTest[i,2] & fusTest[i,3]<1e-2){
                        fusTest[i,4]=3
                }
                if (fusTest[i,1]<fusTest[i,2] & fusTest[i,3]<1e-2){
                        fusTest[i,4]=1
                }
                if (fusTest[i,1]<fusTest[i,2] & fusTest[i,3]>1e-2){
                        fusTest[i,4]=2
                }
                if (fusTest[i,1]>fusTest[i,2] & fusTest[i,3]>1e-2){
                        fusTest[i,4]=2
                }
        }
        print(head(fusTest))

        arrangedfus = matrix(data = NA, nrow = 64, ncol = 2, byrow = FALSE, dimnames = NULL)

        for (i in 1:64) {
                if ( length(which (rownames(fusTest)==codonw[i,]) )<1  ) {
                        arrangedfus[i,1]=codonw[i,]
                        arrangedfus[i,2]=2
                }
                else {
                        arrangedfus[i,1]=rownames(fusTest)[which (rownames(fusTest)==codonw[i,])]
                        arrangedfus[i,2]=as.numeric(fusTest[ which (rownames(fusTest)==codonw[i,]),4])
                }
        }
        fusCOA =matrix(as.numeric(arrangedfus[,2]), nrow = 4, ncol = 16, byrow = TRUE, dimnames = NULL)
        print("Writing results to files!")

        #write.table(data.frame(arrangedfus), file = "fus_arr_new.txt",append = FALSE, sep = ",",row.names = FALSE, col.names = FALSE,qmethod = "escape")
        #write.table(fusCOA, file = "fus_new.coa",append = FALSE, sep = ",",row.names = FALSE, col.names = FALSE)

        return(fusCOA)
}
