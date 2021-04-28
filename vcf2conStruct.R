
direct.vcf2af<-function(input){
  x<-input@gt[,-1]
  x[grep("0/0",x)]<-0
  x[grep("0/1",x)]<-0.5
  x[grep("1/0",x)]<-0.5
  x[grep("1/1",x)]<-1
  x[grep(".:.:.",x)]<-NA
  x<-apply(x,2,as.numeric)
  rownames(x)<-input@fix[,"ID"]
  return(t(x))
}

vcf2conStruct<-function(vcf,metadatafile,outname){
  v<-read.vcfR(vcf)
  meta<-read.table(metadatafile,header=F)
  names(meta)<-c("sample","pop","longitude","latitude")
  v2<-v[,c(1,which(colnames(v@gt)%in%meta$sample))]
  af<-direct.vcf2af(v2)
  i<-data.frame(sample=as.factor(rownames(af)))
  mystrata<-left_join(meta,i,by="sample")
  names(mystrata)
  conStruct_object<-list(allele.frequencies=af,
                         coords=mystrata[,c("longitude","latitude")],
                         geoDist=as.matrix(
                           rdist.earth(mystrata[,c("longitude","latitude")],miles=F)
                         )
  )
  save(conStruct_object,file=paste(outname,".RData",sep=""))
}

require(argparser)

p <- arg_parser("convert VCF to conStruct. Outputs <outname>.RData file. \n
                To load into R for conStruct analysis use load('<outname>.RData').\n
                The object will be named conStruct_object \n\n
                Resuires the following packages: vcfR, adegenet, fields, dplyr, argparser.\n
                Run using Rscript, e.g.\n
                Rscript vcf2conStruct.R --vcf my_file.vcf.gz --metadatafile my_meta.txt --outname my_conStruct",hide.opts = T)

# Add a positional argument
p <- add_argument(p, "--vcf", help="vcf file (can be [b]gzipped)")
p <- add_argument(p, "--metadatafile", help="file with *no header* with columns in this order: 'sample population longitude latitude'")
p <- add_argument(p, "--outname", help="<name of output>.RData")


# parse
args<-parse_args(p)

if(anyNA(c(args$vcf,args$metadatafile,args$outname))) {
  # Print the help message
  print(p)
}else{
  require(vcfR)
  require(fields)
  require(dplyr)
  vcf2conStruct(args$vcf,args$metadatafile,args$outname)
}
















