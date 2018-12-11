
##############################################
# libraries:
##############################################

pkgLoad("lsa") # cosine()
pkgLoad("randomForest") # predict()

options(scipen = 999)


##############################################
# define data files:
##############################################

work_dir <- getwd()
out_dir <- paste0(work_dir, "/results/")
data_dir <- paste0(work_dir, "/data/")
hic_dir <- paste0(work_dir, "/data/dumped/")

hic_domains <- paste0(data_dir, "mouse_limb_E11.5_mm9_mapq30_KR_250kb_all_domains.bed")
pos_set <- paste0(data_dir, "positive_set.bed")
neg_set <- paste0(data_dir, "negative_set.bed")


##############################################
# custom functions:
##############################################

source(paste0(work_dir, "/functions.R"))

##############################################
# some thresholds:
##############################################

binSize <- 10000	#bin sized used to create the matrix
probs.thres <- 0.98	#threshold for upper quantile for hic matrix entries	
cutoff.random.forest <- 0.5 # cutoff for (random forest) predictions
n.tree <- 100 # number of trees in the random forest

##############################################
# Read called domains/TADs:
##############################################

tad <- read.table(file = hic_domains)
colnames(tad) = c("chr","x1", "x2")

##############################################
# fit model (with trainings data):
##############################################

# read hic file:
chr <- "chr1"
hic_file <- list.files(hic_dir, pattern = "chr1_", full.names = T)
m <- read_hic_file(hic_file, binSize, probs.thres)

# get training set:
train_pos <- read.table(file = pos_set)
train_pos <- train_pos[,c(1:3)]
colnames(train_pos) [1:3] <- c("chr","x1", "x2")
train_pos$x1 <- train_pos$x1 - 1
train_pos$class <- 1

train_neg <- read.table(file = neg_set)
train_neg <- train_neg[,c(1:3)]
colnames(train_neg) [1:3] <- c("chr","x1", "x2")
train_neg$x1 <- train_neg$x1 - 1
train_neg$class <- 0

# restrict tads to training chromosome:
this_tad <- tad[which(tad[,1] == chr), c("chr","x1", "x2")]
this_tad$class <- NA

# combine:
df_tad <- rbind(train_pos, train_neg, this_tad)

# define start (b) and end positions (e) in matrix:
df_tad$b <- round((df_tad$x1/binSize + 1),0)
df_tad$e <- round(df_tad$x2/binSize,0)
df_tad <- df_tad[which((df_tad$e - df_tad$b > 2)),]

# initialize metrics:
df_tad$imbalance <- NA
df_tad$imbalance_diff <- NA
df_tad$cosine <- NA
df_tad$direction <- NA

# get metrics:
df_tad_final <- get_sim(df_tad, m, binSize, flex.bins = 5)
df_tad_final <- df_tad_final[which(!is.na(df_tad_final$class)),]

# fit the model:
model <- fit_model(df_tad_final, n.tree = n.tree)


##############################################
# predict imbalanced structures:
##############################################

# get all dumped hic files:
hic_files = list.files(hic_dir, full.names = T)

# predict for each file separately:
result <- c()
for (this_file in hic_files) {
  
  message(paste0("Processing file ",this_file))
  this_chr <- gsub("_.*","",basename(this_file))
  m <- read_hic_file(this_file, binSize, probs.thres)
  
  # restrict tads to training chromosome:
  this_tad <- tad[which(tad[,1] == chr), c("chr","x1", "x2")]
  this_tad$class <- NA
  
  # define start (b) and end positions (e) in matrix:
  this_tad$class <- NA
  this_tad$b <- round((this_tad$x1/binSize + 1),0)
  this_tad$e <- round(this_tad$x2/binSize,0)
  this_tad <- this_tad[which((this_tad$e - this_tad$b > 2)),]
  
  # initialize metrics:
  this_tad$imbalance      <- NA
  this_tad$imbalance_diff <- NA
  this_tad$cosine <- NA
  this_tad$direction <- NA
  
  # column for random forest probabilities:
  this_tad$probability <- NA
  
  # get metrics:
  df_tad_final <- get_sim(this_tad, m, binSize, flex.bins = 5)

  # random forest:
  res <- run_tree(model, df_tad_final)
  idx <- (res[,2] >= cutoff.random.forest)	
  
  if (length(idx) > 0) {
    idx[is.na(idx)] = F
    idx[which(idx == T)] = 1
    idx[which(idx == F)] = 0
    df_tad_final$class = idx
    
    df_tad_final$probability = res[,2]
    df_tad_final[,1] <- paste0("chr",gsub("chr","",df_tad_final[,1]))
    result <- rbind(result, df_tad_final)
  }
}


##############################################
# save results:
##############################################

result <- result[which((as.numeric(result$tuned.x2) - as.numeric(result$tuned.x1)) > 0), ]

out.file <- paste0(out_dir, "Imbalance_quantile.tres",probs.thres,"_random_forest.",cutoff.random.forest,".bed")
write.table(result,
            file = out.file,
            quote = F,
            row.names = F,
            col.names = T,
            sep = "\t")
