# subset the (corrected) CMAP instance file to select only
# relevant entries from the dataset as input for normalization, etc

library(dplyr)
library(data.table)

# input format (delim = ,)
# instance_id,batch_id,cmap_name,INN1,concentration (M),duration (h),cell2,array3,perturbation_scan_id,vehicle_scan_id4,scanner,vehicle,vendor,catalog_number,catalog_name
# output format (delim = \t)
# instance_id	batch_id	cmap_name	INN1	concentration (M)	duration (h)	cell2	array3	scanner	vehicle	vendor	catalog_number	catalog_name	type	scan_id

path_to_instances <- "path/to/instance/file.csv"
#output
path_to_subset_instances <- "path/to/subset/instance/file.tsv"
# control which kind of celllines/microchips you want
celline <- "MCF7"
chip <- "HT_HG-U133A"

instances <- fread(path_to_instances)
long_instance <- melt(instances, measure.vars = c("perturbation_scan_id", "vehicle_scan_id4"),
                      value.name = "scan_id", variable.name = "type") %>%
  unique()

# only subset by cellline & chip
sub_instance <- long_instance[cell2==celline&array3==chip,]
write.table(sub_instance, file = path_to_subset_instances,
            quote=FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# subset for each drug (does not include control smaples!)
path_to_drug_inst <- "path/to/instances/by/drug/folder/"
cmap_names <- instances[cell2==celline&array3==chip,]$cmap_name %>% unique()
cmap_names <- sub("/", "", cmap_names)

lapply(cmap_names, function(n){
  long_instance <- long_instance[cell2==celline&array3==chip&cmap_name==n&type=="perturbation_scan_id",]
  path_to_file <- paste0(path_to_drug_inst, n, "_cels.tsv")
  write.table(long_instance, file=path_to_file,
              sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
  
})