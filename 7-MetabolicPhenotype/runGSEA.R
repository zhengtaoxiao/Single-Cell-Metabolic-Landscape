#########################################################################################################################
runGSEA<-function(expr_data,covariate,test,control,base.name,gmt.file,outdir){
   testname<-paste(base.name,test,'vs',control,sep='_')
   testname<-gsub(' ','_',testname)
   rdata<-expr_data[,colData(expr_data)[[covariate]]%in%c(test,control)]
   rdata<-rdata[,c(grep(paste('^',test,'$',sep=''),colData(rdata)[[covariate]]),
                   grep(paste('^',control,'$',sep=''),colData(rdata)[[covariate]]))]
   
   # save expression matrix
   write.table(rbind(c('symbols',colnames(rdata)),
                     cbind(rownames(rdata),assay(rdata,"exprs"))),
               file='expr.txt',
               quote=F,
               sep='\t',
               col.names=F,
               row.names=F)
   
   #save the cls file
   pheno<-as.character(colData(rdata)[[covariate]])
   con<-file('pheno.cls',open='w')
   write(paste(length(pheno),'2 1'),con)
   write(paste('# ',test,' ',control,sep=''),con)
   classes<-''
   for (i in 1:length(pheno)){
      classes<-paste(classes,pheno[i])
   }
   write(classes,con)
   close(con)         
   
   #call java gsea version 
   command <- paste('java -Xmx512m -cp ../gsea-3.0.jar xtools.gsea.Gsea -res expr.txt -cls pheno.cls#',test,'_versus_',control,' -gmx ',gmt.file,
                ' -collapse false -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label ',testname, 
                ' -metric Diff_of_Classes -sort real -order descending -include_only_symbols false -make_sets true -median false -num 100',
                ' -plot_top_x 20 -rnd_seed 123456 -save_rnd_lists false -set_max 10000 -set_min 5 -zip_report false -out ', outdir, ' -gui false',sep='')
   if(get_os() == "win"){
      system(command,show.output.on.console=F)
   }else{
      system(command)
   }
   
   unlink(c('pheno.cls','expr.txt'))
}
