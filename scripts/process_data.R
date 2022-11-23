outdir = 'results/input_data/'
dir.create(outdir, showWarnings = F)

### process cytof data
print('### process cytof data')
cytof_fn = 'data/cmi_pb/cytof.2020.pma.day0.proc.tsv'
cytof_new = paste0(outdir, 'cytof.2020.momix.day0.input.tsv')
cytof = read.table(cytof_fn, sep='\t')
cytof = t(cytof)
#cytof = na.omit(cytof)
cytof = cytof[rowSums(is.na(cytof)) == 0, colSums(is.na(cytof)) == 0]
cytof = cytof[rowSums(cytof == "") == 0, colSums(cytof== "") == 0]

cytof[1,1] = "probe"
write.table(cytof,
            file=cytof_new, 
            quote = F,
            row.names = F,
            sep = '\t',
            col.names = F)

print('Sanity check -- cytof has the following dims:')
str(dim(cytof))

### process olink data
print('### process olink data')
olink_fn = 'data/cmi_pb/olink.2020.pma.day0.proc.tsv'
olink_new = paste0(outdir, 'olink.2020.momix.day0.input.tsv')
olink = read.table(olink_fn, sep='\t')
olink = t(olink)
olink[1,1] = "probe"
#olink = olink[rowSums(is.na(olink)) == 0, colSums(is.na(olink)) == 0]
print('before')
print(olink)

olink = olink[, colSums(olink == "") < 5]
print('after')
print(olink)


olink = olink[rowSums(olink == "") == 0,]
print('after2')
print(olink)

write.table(olink,
            file=olink_new, 
            quote = F,
            row.names = F,
            sep = '\t',
            col.names = F)
print('Sanity check -- olink has the following dims:')
str(dim(olink))

### process rnaseq data
print('### process rnaseq data')
rnaseq_fn = 'data/cmi_pb/rnaseq.2020.pma.day0.proc.tsv'
rnaseq_new = paste0(outdir, 'rnaseq.2020.momix.day0.input.tsv')
rnaseq = read.table(rnaseq_fn, sep='\t')
rnaseq = t(rnaseq)
rnaseq[1,1] = "probe"
rnaseq = rnaseq[rowSums(is.na(rnaseq)) == 0, colSums(is.na(rnaseq)) == 0]
#rnaseq = na.omit(rnaseq)
write.table(rnaseq,
            file=rnaseq_new, 
            quote = F,
            row.names = F,
            sep = '\t',
            col.names = F)

print('Sanity check -- rnaseq has the following dims:')
str(dim(rnaseq))

