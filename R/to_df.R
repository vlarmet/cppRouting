to_df<-function(Graph){
  dat<-Graph$data
  dat$from<-Graph$dict$ref[match(dat$from,Graph$dict$id)]
  dat$to<-Graph$dict$ref[match(dat$to,Graph$dict$id)]
  return(dat)
}
