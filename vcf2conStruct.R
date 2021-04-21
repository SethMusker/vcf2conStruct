

# functions
vcf2gi<-function(vcf,metadatafile){
  ## metadatafile must be 'sample pop long lat'
  g<-vcfR2genind(read.vcfR(vcf))
  meta<-read.table(metadatafile,header=F)
  names(meta)<-c("sample","pop","longitude","latitide")
  g2<-g[indNames(g)%in%meta$sample,]
  i<-data.frame(sample=indNames(g2))
  mystrata<-left_join(i,meta,by="sample")
  pop(g2)<-mystrata$pop
  strata(g2)<-mystrata
  g2@other$latlong<-mystrata[,c("longitude","latitide")]
  return(g2)
}

genind2conStruct<-function(genind){
  df<-as.matrix(genind2df(genind))
  df[is.na(df)]<-"-9"
  df[df=="00"]<-0
  df[df=="01"|df=="10"]<-0.5
  df[df=="11"]<-1
  df[df=="-9"]<-NA
  conStruct.manual<-as.data.frame(df)
  conStruct.data<-list(allele.frequencies=apply(as.matrix(conStruct.manual[,-1]),2,as.numeric))
  rownames(conStruct.data$allele.frequencies)<-rownames(df)
  conStruct.data$coords<-as.matrix(genind@other$latlong)
  conStruct.data$geoDist<-as.matrix(rdist.earth(conStruct.data$coords,miles=F))
  return(conStruct.data)
}

vcf2conStruct<-function(vcf,metadatafile,outname){
  conStruct_object<-genind2conStruct(vcf2gi(vcf,metadatafile))
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
  require(adegenet)
  require(fields)
  require(dplyr)
  vcf2conStruct(args$vcf,args$metadatafile,args$outname)
}
















