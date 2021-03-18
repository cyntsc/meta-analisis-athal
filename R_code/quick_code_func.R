## Quick functions

library(here);   
here();
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

## Open a PNG object
library(png)
img1_path <- "/data/run/cyntsc/Project_athal_wgcna/results/plots/Athal_GO.png"
img1 <- readPNG(img1_path, native = TRUE, info = TRUE)
attr(img1, "info")

## Clear all variables created }
rm(list = ls());


## retrieve column names
names(dataf)
## retrieve values of a column
dataf$column_name
## Sets a column name as () index row.names )
athalData <- read.csv(here("data", "all_log2_tidy.csv"), header=TRUE, row.names='Genes', sep='\t')

# Convert a df in matriz numeric
is.matrix(datExpr);
class(athalData);
athalData = data.matrix(athalData, rownames.force = NA)
is.matrix(athalData);
