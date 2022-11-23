outdir = 'results/input_data/'
dir.create(outdir, showWarnings = F)

### process cytof data
print('### process cytof data')
cytof_fn = 'data/cmi_pb/cytof.pma.day0.proc.tsv'
cytof_new = paste0(outdir, 'cytof.momix.day0.input.tsv')
cytof = read.table(cytof_fn)
cytof = t(cytof)
cytof[1,1] = "probe"
write.table(cytof,
            file=cytof_new, 
            quote = F,
            row.names = F,
            col.names = F)

print('Sanity check -- cytof has the following dims:')
str(dim(cytof))

### process olink data
print('### process olink data')
olink_fn = 'data/cmi_pb/olink.pma.day0.proc.tsv'
olink_new = paste0(outdir, 'olink.momix.day0.input.tsv')
olink = read.table(olink_fn)
olink = t(olink)
olink[1,1] = "probe"
write.table(olink,
            file=olink_new, 
            quote = F,
            row.names = F,
            col.names = F)
print('Sanity check -- olink has the following dims:')
str(dim(olink))

### process rnaseq data
print('### process rnaseq data')
rnaseq_fn = 'data/cmi_pb/rnaseq.pma.day0.proc.tsv'
rnaseq_new = paste0(outdir, 'rnaseq.momix.day0.input.tsv')
rnaseq = read.table(rnaseq_fn)
rnaseq = t(rnaseq)
rnaseq[1,1] = "probe"
write.table(rnaseq,
            file=rnaseq_new, 
            quote = F,
            row.names = F,
            col.names = F)

print('Sanity check -- rnaseq has the following dims:')
str(dim(rnaseq))

