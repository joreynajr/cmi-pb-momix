## Support function to apply log2(+1) to a matrix
log2matrix <- function(folder, file.name){
    
    # Read table
    fn <- paste(folder, file.name, sep="/")
    data <- as.matrix(read.table(fn, sep=" ", row.names=1, header=TRUE))
    
    # Apply transformation
    data <- log2(data + 1)
    
    # Output file name
    output <- paste(folder, paste0("log_",file.name), sep="/")
    
    # Export transformed data
    write.table(data, output, sep=" ", col.names=TRUE, row.names=TRUE)
    
    # ?
    # Adding probe for some reason? Probably just tries to ensure
    # that the index column has a name which is not saved by default in R.
    # Modification: sed works differently on MAC's then on Linux.
    # Previously this function used sed -i which means inplace for the Linux 
    # version but serves as some extension modifier in the Mac version.
    tmp = paste0(output, ".tmp")
    cmd = paste("sed '1s/^/probe\t/'", output, ">", tmp, sep=" ")
    print(cmd)
    system(cmd)
    file.rename(tmp, output)
}