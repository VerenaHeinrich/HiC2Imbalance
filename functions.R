

#############################
# check packages:
#############################

pkgTest <- function(pkg){
  if (!pkg %in% installed.packages()) {
    stop(paste0("package ",pkg, " not found."),call. = FALSE)
  }
} 

#############################
# load packages:
#############################

pkgLoad <- function(pkg) {
  pkgTest(pkg)
  cat(paste0(".. load package ", pkg))
  suppressMessages(library(pkg, character.only = TRUE))
}

#############################
# read hic file:
#############################

read_hic_file <- function(hic_file, binSize, probs.thres){
  
  tab <- read.table(file = gzfile(hic_file), col.names = c("x", "y", "count"))
  tab$count[which(is.na(tab$count))] <- 0	
  
  tab$i <- tab$x/binSize + 1
  tab$j <- tab$y/binSize + 1
  
  thresh <- quantile(x = tab$count, probs = probs.thres, na.rm = T)
  tab$count[tab$count > thresh] <- thresh 
  
  n <- max(c(tab$i, tab$j))
  m <- matrix(0, ncol = n, nrow = n)
  
  # Make the matrix symmetic
  m[cbind(tab$i, tab$j)] <- tab$count
  m[cbind(tab$j, tab$i)] <- tab$count
  
  # remove all values < 0
  m[m < 0] <- 0  
  
  return(m)
}

#############################
# get inbalance metrics:
#############################

get_sim <- function(df_tad, m, binSize, flex.bins){

	#define fine tuned tads:
	tads.tuned = c()

	for (i in 1:(dim(df_tad)[1])) {

		# define start end end coordinates
		this.start <- df_tad$b[i] - flex.bins
 		this.end <- df_tad$e[i] + flex.bins

		#if subscript out of bounds:
		if (this.start < 1) {
			this.start <- df_tad$b[i] - df_tad$b[i] - 1
		}
	  if (this.end > dim(m)[1]) {
			this.end <- df_tad$e[i] + (dim(m)[1] - df_tad$e[i])
		}
		if (this.end <= this.start) {
			this.start <- df_tad$b[i]
			this.end <- df_tad$e[i]

			new_start <- this.start
			new_end <- this.end

			this.m <- m[this.start:this.end,this.start:this.end]
			dim <- dim(this.m)[2]

			this.v1.max.diff <- this.m[1,]
			this.v2.max.diff <- rev(this.m[,dim])
		}else{

			#define total number of bins:
			flex.bins.final <- flex.bins*2 + 1

			# define sub-matrix:
			this.m <- m[this.start:this.end,this.start:this.end]
			dim <- dim(this.m)[2]

			# define vectors at each side of the upper triangle:
			this.v1 <- this.m[1:flex.bins.final,1:(dim - 1)]
			this.v2 <- this.m[2:dim,(dim - flex.bins.final + 1):dim]

			# get vector with highest mean:
			this.v1.max.diff.idx <- order(apply(this.v1, 1, function(x) max(mean(x,na.rm = T),na.rm = T)),decreasing = T)[1]
			this.v2.max.diff.idx <- order(apply(this.v2, 2, function(x) max(mean(x,na.rm = T),na.rm = T)),decreasing = T)[1]

			max.na <- max(this.v1.max.diff.idx,(1 - this.v1.max.diff.idx))
			this.v1.max.diff <- this.v1[this.v1.max.diff.idx,max.na:dim(this.v1)[2]]
			this.v2.max.diff <- rev(this.v2[,this.v2.max.diff.idx])[max.na:dim(this.v2)[1]]

			# get fine tuned domains;
			new_start <- this.start + this.v1.max.diff.idx - 1
			new_end <- this.end - flex.bins.final + this.v2.max.diff.idx - 1
		}

		tads.tuned <- rbind(  tads.tuned, 
                	        c(as.character(df_tad[1,1]),
                	          new_start*binSize,
                	          new_end*binSize))


 		#define metrics:
	  df_tad$imbalance[i] <- abs(log2(mean(this.v1.max.diff)/mean(this.v2.max.diff)))
		df_tad$imbalance_diff[i] <- abs(mean(this.v1.max.diff) - mean(this.v2.max.diff))
  	df_tad$direction[i] <- (log2(mean(this.v1.max.diff)/mean(this.v2.max.diff))) > 0

		#cosine based on similarity on matrix:
  	df_tad$cosine[i] <- abs(cosine(this.v1.max.diff,this.v2.max.diff))
	}

	tads.tuned[,2] <- gsub(" ","",format(as.numeric(tads.tuned[,2]),scientific = F))
	tads.tuned[,3] <- gsub(" ","",format(as.numeric(tads.tuned[,3]),scientific = F))
	df_tad$tuned.x1 <- tads.tuned[,2]
	df_tad$tuned.x2 <- tads.tuned[,3]

	return(df_tad)
}

#############################
# fit the model:
#############################

fit_model <- function(data, n.tree){
  
  model <- randomForest(as.factor(class) ~ imbalance + imbalance_diff  + cosine,
                     data = data,
                     ntree = n.tree)
  return(model)
}


#############################
# get predicted values:
#############################

run_tree <- function(model, df_tad_final){
	
	res <- predict(model,
	               df_tad_final,
	               type = "prob")
	return(res)
}
