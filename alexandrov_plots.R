# If middle base is A or G, change to c or t and take reverse complement

# FIRST
# maftool get_context maf [num_bases] > new_maf

# Then take output maf:

# count num muts occur in each trinucleotide in sample
# ---------------
# how many total occurences of the trinucleotide in exome

library(dplyr)


trim = function (x){ gsub("^\\s+|\\s+$", "", x) }

base_complement = function(base){
  if(base == 'A'){ return('T')}
  else if(base == 'C'){ return('G')}
  else if(base == 'T'){ return('A')}
  else if(base == 'G'){ return('C')}
  else{ stop(paste('Not a base:', base))}
}

does_this_need_to_be_reversed = function(base){
  if(base == 'A' | base == 'G'){ return('yes') }
  else{ return('no') }
}

folded_tri_context = function( tris ){
  folded = unlist(sapply(tris, function(x){
    tri_split = strsplit(x, '')[[1]]
    if(does_this_need_to_be_reversed(tri_split[2]) == 'yes'){
      tri_split_rev_comp = paste( base_complement(tri_split[3]), 
                                  base_complement(tri_split[2]), 
                                  base_complement(tri_split[1]), 
                                  sep='')
      return(tri_split_rev_comp)
    }
    else{ return(x) }
  }))
  return(folded)
}

group_exome_counts = function(ec){
  library(dplyr)
  ec = group_by(ec, halved_context)
  x = dplyr:::summarise(ec, exome_trinuc_count = sum(count))
  return(x)
}

make_full_context = function(x_context){
  x = apply(x_context, 1, function(a){
    splt = strsplit(a['context'], '')[[1]]
    return(paste( splt[1], a['Reference_Allele'], splt[3], sep=''))
    })
  return(x)
}

make_bracket_columns = function(maf_cols){
  idk = apply(maf_cols, 1, function(x){
    split_context = strsplit(x['full_context'], '')[[1]]
    if( x['Reference_Allele'] == 'A' | x['Reference_Allele'] == 'G'){
      return( paste( base_complement(split_context[3]), '[', base_complement(x['Reference_Allele']), '/', base_complement(x['Tumor_Seq_Allele1']), ']', base_complement(split_context[1]) , sep = ''))
    }
    else{
      return( paste( split_context[1], '[',x['Reference_Allele'], '/', x['Tumor_Seq_Allele1'], ']', split_context[3], sep='' ))
    }
  })
  return(idk)
}

add_context_columns = function(maf_loc, scratch_maf_loc){
  system(paste("maftool get_context", maf_loc, "1 >", scratch_maf_loc, sep = ' '))
  maf = read.table(scratch_maf_loc, sep='\t', header = T, stringsAsFactors = F)
  maf = maf[ which(maf$Variant_Type == 'SNP'), ]

  maf$full_context = make_full_context(maf[, c('Reference_Allele', 'context')])
  maf$halved_context = folded_tri_context(maf$full_context)
  maf$tri_context = make_bracket_columns(maf[ , c('Reference_Allele', 'Tumor_Seq_Allele1', 'full_context', 'halved_context')])
  return(maf)
}

graph_setup = function(count_context_table){

}

do_it = function(maf_loc, scratch_maf_loc){
  #filter out non-SNPs
  maf = add_context_columns(maf_loc, scratch_maf_loc)
  
  sample_trinuc_counts = as.data.frame(table(maf$halved_context))
  colnames(sample_trinuc_counts) = c('halved_context', 'sample_trinuc_count')

  exome_count = read.table('/home/chapmanjs/tricount/exon_triCount.txt', sep='\t', header = T, stringsAsFactors = F)
  exome_count$halved_context = folded_tri_context(exome_count$tri)
  exome_count_folded = group_exome_counts(exome_count)
  new_table = as.data.frame(table(maf$tri_context), stringsAsFactors = F)
  colnames(new_table) = c('tri_context', 'count')
  new_table$tri_base = unlist(sapply(new_table$tri_context, function(x){
    splt = strsplit(x, '')[[1]]
    return(paste(splt[c(1,3,7)], collapse = ''))
    }))
  new_table$exome_count = sapply(new_table$tri_base, function(x){
    ind = match(x, exome_count_folded$halved_context)
    return(exome_count_folded[ind, 'exome_trinuc_count'])
    })

  new_table$ratio = new_table$count / new_table$exome_count
  mut_count_sum = sum(new_table$count)
  new_table$normalized_ratio = new_table$count / mut_count_sum

  all_96 = read.table('~/all_contexts.txt', sep='\n', header = F, stringsAsFactors = F)
  all_96 = all_96$V1
 
  missing_trinuc_contexts = all_96[(!(all_96 %in% new_table$tri_context))]
  for( i in missing_trinuc_contexts){
    new_table = rbind(new_table, c(i, 0, 'XXX', 0, 0, 0))
  }

  library(plyr)
  # new_table$tri_context <- as.factor(new_table$tri_context)
  new_table$color_ctxt = unlist(sapply(new_table$tri_context, function(x){
    splt = strsplit(x, '')[[1]]
    return(paste(splt[3], '>', splt[5], sep=''))
    }))
  new_table = new_table[ order(new_table$color_ctxt), ]

  new_table$tri_context  = factor(new_table$tri_context, ordered = T)
  
  c(rep(1, 16), rep(2, 16), rep(3, 16), rep(4, 16), rep(5, 16), rep(6, 16))
  # new_table$color_ctxt <- revalue(new_table$tri_context, c("A[C/A]A"="C>A", "A[C/A]C"="C>A", "A[C/A]G"="C>A", "A[C/A]T"="C>A", "C[C/A]A"="C>A", "C[C/A]C"="C>A", "C[C/A]G"="C>A", "C[C/A]T"="C>A",
  #                                                         "G[C/A]A"="C>A", "G[C/A]C"="C>A", "G[C/A]G"="C>A", "G[C/A]T"="C>A", "T[C/A]A"="C>A", "T[C/A]C"="C>A", "T[C/A]G"="C>A", "T[C/A]T"="C>A", "A[C/G]A"="C>G", "A[C/G]C"="C>G", "A[C/G]G"="C>G",
  #                                                         "A[C/G]T"="C>G", "C[C/G]A"="C>G", "C[C/G]C"="C>G", "C[C/G]G"="C>G", "C[C/G]T"="C>G", "G[C/G]A"="C>G", "G[C/G]C"="C>G", "G[C/G]G"="C>G", "G[C/G]T"="C>G", "T[C/G]A"="C>G", "T[C/G]C"="C>G",
  #                                                         "T[C/G]G"="C>G", "T[C/G]T"="C>G", "A[C/T]A"="C>T", "A[C/T]C"="C>T", "A[C/T]G"="C>T", "A[C/T]T"="C>T", "C[C/T]A"="C>T", "C[C/T]C"="C>T", "C[C/T]G"="C>T", "C[C/T]T"="C>T", "G[C/T]A"="C>T",
  #                                                         "G[C/T]C"="C>T", "G[C/T]G"="C>T", "G[C/T]T"="C>T", "T[C/T]A"="C>T", "T[C/T]C"="C>T", "T[C/T]G"="C>T", "T[C/T]T"="C>T", "A[T/A]A"="T>A", "A[T/A]C"="T>A", "A[T/A]G"="T>A", "A[T/A]T"="T>A",
  #                                                         "C[T/A]A"="T>A", "C[T/A]C"="T>A", "C[T/A]G"="T>A", "C[T/A]T"="T>A", "G[T/A]A"="T>A", "G[T/A]C"="T>A", "G[T/A]G"="T>A", "G[T/A]T"="T>A", "T[T/A]A"="T>A", "T[T/A]C"="T>A", "T[T/A]G"="T>A",
  #                                                         "T[T/A]T"="T>A", "A[T/C]A"="T>C", "A[T/C]C"="T>C", "A[T/C]G"="T>C", "A[T/C]T"="T>C", "C[T/C]A"="T>C", "C[T/C]C"="T>C", "C[T/C]G"="T>C", "C[T/C]T"="T>C", "G[T/C]A"="T>C", "G[T/C]C"="T>C",
  #                                                         "G[T/C]G"="T>C", "G[T/C]T"="T>C", "T[T/C]A"="T>C", "T[T/C]C"="T>C", "T[T/C]G"="T>C", "T[T/C]T"="T>C", "A[T/G]A"="T>G", "A[T/G]C"="T>G", "A[T/G]G"="T>G", "A[T/G]T"="T>G", "C[T/G]A"="T>G",
  #                                                         "C[T/G]C"="T>G", "C[T/G]G"="T>G", "C[T/G]T"="T>G", "G[T/G]A"="T>G", "G[T/G]C"="T>G", "G[T/G]G"="T>G", "G[T/G]T"="T>G", "T[T/G]A"="T>G", "T[T/G]C"="T>G", "T[T/G]G"="T>G", "T[T/G]T"="T>G"))

  library(ggplot2)
  
  ggp <- ggplot(new_table, aes(x=tri_context, y=normalized_ratio, fill=color_ctxt))
  print(ggp + geom_histogram(stat="identity", drop=F) + scale_fill_manual(values=c("#1EBBEB", "#000000", "#E0242A", "#A3A3A3", "#A3CC6E", "#EDC6C2")) + theme_bw() + theme(axis.text.x=element_text(angle = 45, vjust = 0.5)))
                            

  return(new_table)
}








