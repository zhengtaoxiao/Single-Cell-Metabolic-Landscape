runGSEA_preRank<-function(preRank.matrix,gmt.file,outname){
   #descending numerical order
   #dump preRank into a tab-delimited txt file
   write.table(preRank.matrix,
               file='prerank.rnk',
               quote=F,
               sep='\t',
               col.names=F,
               row.names=T)
      
   #call java gsea version
   system(paste('java -Xmx512m -cp ../gsea-3.0.jar xtools.gsea.GseaPreranked -gmx ', gmt.file, ' -norm meandiv -nperm 1000 -rnk prerank.rnk ',
                ' -scoring_scheme weighted -make_sets true -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false ',
                ' -out preRankResults -create_svgs true -gui false -rpt_label ',outname, sep=''), show.output.on.console=F)
  unlink(c('prerank.txt'))
}
